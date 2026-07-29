/* fairchild — C API for the fairchild analog + photonic circuit simulator.
 *
 * Two layers over the same solver:
 *
 *   Batch     fc_sim_new -> fc_load_* -> fc_run_tran -> fc_signal
 *   Stepping  fc_stepper_new -> { fc_get_node / fc_set_source / fc_step }*
 *
 * The stepping layer holds a transient open between timesteps so the host
 * program drives the clock — that is the mixed-signal path, with signals going
 * both ways every step.
 *
 * Link against libfairchild_c (shared or static):
 *     cargo build -p fairchild-c --release
 *     cc prog.c -I crates/fairchild-c/include -L target/release -lfairchild_c
 *
 * The "_c" suffix keeps the artifact from colliding with the Python extension
 * module, which also builds as libfairchild.  Rename or symlink it freely.
 *
 * Conventions:
 *   - Functions returning int return FC_OK or an FC_ERR_* code; fc_error /
 *     fc_stepper_error give the message, owned by the handle and valid until
 *     the next call on it.
 *   - NULL handles and NULL string arguments return FC_ERR_ARG, never crash.
 *   - A handle is not thread-safe.  Use one per thread; there is no global
 *     state, so independent simulations run concurrently without locking.
 */
#ifndef FAIRCHILD_H
#define FAIRCHILD_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FC_OK            0
#define FC_ERR_ARG       1  /* NULL pointer, bad UTF-8, or no netlist loaded  */
#define FC_ERR_PARSE     2  /* netlist could not be read or parsed            */
#define FC_ERR_SIM       3  /* no convergence, singular matrix, floating node */
#define FC_ERR_NOT_FOUND 4  /* no such node, source, or signal                */
#define FC_ERR_PANIC     5  /* internal fault: free the handle, do not reuse  */

typedef struct FcSim     fc_sim;
typedef struct FcStepper fc_stepper;

/* Library version string; static, never NULL. */
const char *fc_version(void);

/* ---- handle lifecycle ------------------------------------------------- */

fc_sim *fc_sim_new(void);
void    fc_sim_free(fc_sim *sim);              /* NULL is a no-op */
/* Last error on this handle, or NULL if the last call succeeded. */
const char *fc_error(const fc_sim *sim);

/* ---- loading and editing the netlist ---------------------------------- */

int fc_load_file(fc_sim *sim, const char *path);
int fc_load_string(fc_sim *sim, const char *text);

/* Solver option, as ".options KEY=VALUE" would set it: reltol, abstol, vntol,
 * gmin, itl1, itl4, maxstep, method (be|tr|gear), solver (dense|sparse|klu),
 * variable_step, uic, temp, lambda_center_nm, ... */
int fc_set_option(fc_sim *sim, const char *key, const char *value);

/* Retarget an element parameter: fc_set_param(sim, "R1", "value", 2e3).
 * Passives take "value" or their physical name; V/I sources take "dc"/"value";
 * MOSFET and OSDI instances take any instance parameter. */
int fc_set_param(fc_sim *sim, const char *element, const char *param, double value);

/* Replace a source waveform with a PWL table of n points (t ascending).
 * Arrays are copied.  Use this when the stimulus is known up front; for a value
 * decided during the run, use fc_set_source on a stepper. */
int fc_set_source_pwl(fc_sim *sim, const char *name,
                      const double *t, const double *v, size_t n);

/* ---- batch analyses --------------------------------------------------- */

int fc_run_op(fc_sim *sim);
int fc_run_tran(fc_sim *sim, double step, double stop);

/* Operating-point values, after fc_run_op. */
int fc_op_node(fc_sim *sim, const char *node, double *value);
int fc_op_current(fc_sim *sim, const char *vsrc, double *value);

/* Transient waveforms, after fc_run_tran.  Name is "time", "V(node)", or
 * "I(vsrc)".  *data borrows the handle's storage — do not free it; it is
 * invalidated by the next run or by fc_sim_free. */
int fc_signal(fc_sim *sim, const char *name, const double **data, size_t *len);

/* Enumerate available signal names.  fc_signal_name returns NULL out of range;
 * the pointer is valid until the next run. */
size_t      fc_signal_count(const fc_sim *sim);
const char *fc_signal_name(const fc_sim *sim, size_t i);

/* ---- host-driven stepping (mixed signal) ------------------------------ */

/* Open a transient with fixed timestep `step`, solving the operating point so
 * the handle starts at t = 0.  Snapshots the netlist: later fc_set_param calls
 * on `sim` do not affect it, and it stays valid after fc_sim_free.
 * Returns NULL on failure — fc_error(sim) says why. */
fc_stepper *fc_stepper_new(fc_sim *sim, double step);
void        fc_stepper_free(fc_stepper *st);   /* NULL is a no-op */
const char *fc_stepper_error(const fc_stepper *st);

/* Advance one timestep.  t_out may be NULL.  On FC_ERR_SIM the stepper still
 * holds the last accepted timepoint, so a host can back off and retry. */
int fc_step(fc_stepper *st, double *t_out);

/* Step until the time reaches t_target.  The step size is fixed, so this lands
 * on the first grid point at or past t_target; t_out (may be NULL) reports
 * where.  Already past t_target => no steps taken. */
int fc_advance_to(fc_stepper *st, double t_target, double *t_out);

double fc_time(const fc_stepper *st);       /* seconds; negative if st is NULL */
double fc_step_size(const fc_stepper *st);

/* Analog -> digital: read the present timepoint. */
int fc_get_node(fc_stepper *st, const char *node, double *value);
int fc_get_current(fc_stepper *st, const char *vsrc, double *value);

/* Digital -> analog: hold a source at `value` from the next step on.  Zero-order
 * hold, matching ngspice's GetVSRCData semantics. */
int fc_set_source(fc_stepper *st, const char *name, double value);

/* Enumerate node names.  fc_node_name copies node i into buf (NUL-terminated,
 * truncated to cap) and returns the length needed excluding the NUL, so passing
 * a NULL buf measures.  Returns -1 if i is out of range. */
size_t    fc_node_count(const fc_stepper *st);
ptrdiff_t fc_node_name(const fc_stepper *st, size_t i, char *buf, size_t cap);

#ifdef __cplusplus
}
#endif

#endif /* FAIRCHILD_H */
