# legacy/

Pre-Phase-B Verilog-A model library and associated tooling.

## va-models/

Verilog-A source files for photonic and electronic device models that were
the original OSDI-based simulation path. After Phase B, all of these models
have native Rust equivalents in `crates/fairchild-core/src/models/`.

**These files are no longer maintained.** They are kept for reference and for
the OSDI loader tests in `crates/fairchild-osdi/tests/`.

### Building the .osdi binaries

The compiled `.osdi` artifacts are not tracked in git. To rebuild them:

```sh
cd legacy/va-models
bash build.sh          # requires openvaf-r on PATH
```

Or individually:

```sh
openvaf-r legacy/va-models/photonic/waveguide.va \
    --output legacy/va-models/build/waveguide.osdi
```
