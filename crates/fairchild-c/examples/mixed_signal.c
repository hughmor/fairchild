/* Mixed-signal co-simulation from C: a digital state machine driving an analog
 * RC load, deciding what to do next from the analog result.
 *
 * The "digital" side here is a clocked comparator + counter — stand-in for
 * whatever event-driven simulator you are integrating.  The pattern is the same
 * regardless: advance the analog side to the next digital clock edge, sample it,
 * decide, drive back.
 *
 *   cargo build -p fairchild-c --release
 *   cc -O2 -o mixed_signal crates/fairchild-c/examples/mixed_signal.c \
 *      -I crates/fairchild-c/include -L target/release -lfairchild_c
 *   ./mixed_signal
 */
#include <stdio.h>
#include <stdlib.h>
#include "fairchild.h"

/* R-C load driven by VDRIVE, sensed at OUT.  tau = 1k * 10p = 10 ns. */
static const char *NETLIST =
    "* digitally driven RC load\n"
    "VDRIVE drv 0 DC 0\n"
    "R1 drv out 1k\n"
    "C1 out 0 10p\n"
    ".tran 1n 1u\n"
    ".end\n";

#define CLOCK_PERIOD 2e-9   /* digital clock: 500 MHz          */
#define ANALOG_STEP  1e-10  /* 20 analog steps per clock edge  */
#define VDD          1.8
#define VTHRESH      0.9
#define N_CYCLES     400

static int die(const char *what, const char *msg) {
    fprintf(stderr, "%s: %s\n", what, msg ? msg : "(no message)");
    return 1;
}

int main(void) {
    printf("fairchild %s\n", fc_version());

    fc_sim *sim = fc_sim_new();
    if (!sim) return die("fc_sim_new", "out of memory");

    if (fc_load_string(sim, NETLIST) != FC_OK)
        return die("fc_load_string", fc_error(sim));

    /* Open the transient.  The operating point is solved here, so the first
     * sample below is the real t=0 state, not a guess. */
    fc_stepper *st = fc_stepper_new(sim, ANALOG_STEP);
    if (!st) return die("fc_stepper_new", fc_error(sim));
    fc_sim_free(sim);  /* the stepper snapshotted the netlist; sim is done */

    /* Digital state: drive high until OUT crosses VTHRESH, then low, and count
     * the crossings.  A real integration would hand this to the digital kernel. */
    int driving_high = 1;
    long transitions = 0;
    double v_out = 0.0;
    /* Track the ripple band: the sampled value alone aliases with the clock. */
    double v_min = 1e30, v_max = -1e30;

    if (fc_set_source(st, "VDRIVE", VDD) != FC_OK)
        return die("fc_set_source", fc_stepper_error(st));

    for (long cycle = 0; cycle < N_CYCLES; cycle++) {
        double t_edge = (cycle + 1) * CLOCK_PERIOD, t_now = 0.0;

        /* --- analog advances to the next clock edge --- */
        if (fc_advance_to(st, t_edge, &t_now) != FC_OK)
            return die("fc_advance_to", fc_stepper_error(st));

        /* --- analog -> digital: sample --- */
        if (fc_get_node(st, "out", &v_out) != FC_OK)
            return die("fc_get_node", fc_stepper_error(st));

        /* --- digital decides --- */
        int want_high = (v_out < VTHRESH);
        if (cycle > 20) { /* skip the initial charge-up ramp */
            if (v_out < v_min) v_min = v_out;
            if (v_out > v_max) v_max = v_out;
        }

        /* --- digital -> analog: drive --- */
        if (want_high != driving_high) {
            driving_high = want_high;
            transitions++;
            if (fc_set_source(st, "VDRIVE", driving_high ? VDD : 0.0) != FC_OK)
                return die("fc_set_source", fc_stepper_error(st));
        }

        if (cycle % 50 == 0)
            printf("  t = %8.2f ns   out = %6.4f V   drive = %s\n",
                   t_now * 1e9, v_out, driving_high ? "HIGH" : "LOW");
    }

    printf("\n%ld clock cycles, %ld drive transitions.\n"
           "out ripples in [%.4f, %.4f] V around the %.2f V threshold "
           "(one clock of RC slew per decision).\n",
           (long)N_CYCLES, transitions, v_min, v_max, VTHRESH);

    if (transitions < 2) {
        fprintf(stderr, "FAIL: the loop never toggled — is the netlist driving?\n");
        fc_stepper_free(st);
        return 1;
    }
    printf("OK: closed-loop bang-bang held the node at its threshold.\n");

    fc_stepper_free(st);
    return 0;
}
