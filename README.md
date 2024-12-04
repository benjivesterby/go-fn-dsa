# FN-DSA (in Go)

FN-DSA is a new *upcoming* post-quantum signature scheme, currently
being defined by NIST as part of their [Post-Quantum Cryptography
Standardization](https://csrc.nist.gov/pqc-standardization) project.
FN-DSA is based on the [Falcon](https://falcon-sign.info/) scheme.

**WARNING:** As this file is being written, no FN-DSA draft has been
published yet, and therefore what is implemented here is *not* the
"real" FN-DSA; such a thing does not exist yet. When FN-DSA gets
published (presumably as a draft first, but ultimately as a "final"
standard), this implementation will be adjusted accordingly.
Correspondingly, it is expected that **backward compatiblity will NOT be
maintained**, i.e. that keys and signatures obtained with this code may
cease to be accepted by ulterior versions. Only version 1.0 will provide
such stability, and it will be published only after publication of the
final FN-DSA standard.

This implementation is the Go variant of the
[C](https://github.com/pornin/c-fn-dsa/) and
[Rust](https://github.com/pornin/rust-fn-dsa/) implementations. All
three implementations are interoperable and mostly feature parity.
This pure Go implementation tends to be about 1.7x slower than the C
and Rust versions, in part because Go does not support SSE2 or NEON
intrinsics. The performance achieved with this code is still quite
acceptable (e.g. it remains substantially faster than RSA-2048 on
the architectures for which hardware floating-point can be safely
leveraged, i.e. 386+SSE2, amd64, arm64 and riscv64).

The API is documented in [doc.go](fndsa/doc.go); top-level functions for
key pair generation, signing, and verifying, are found in
[kgen.go](fndsa/kgen.go), [sign.go](fndsa/sign.go), and
[vrfy.go](fndsa/vrfy.go), respectively. A few extra public definitions
are in [util.go](fndsa/util.go).

Internally, the three components of the scheme are segregated, and
(for instance) a reduced implementation that supports only signature
verification can be obtained by simply deleting the files that are
relevant to key pair generation and to signature generation:

  - All files `kgen*.go` are used only for key pair generation.

  - All files `sign*.go` are used only for signature generation.

  - All files `vrfy*.go` are used only for signature verification.

  - [fndsa_test.go](fndsa/fndsa_test.go) is a test file that exercises
    the full scheme, and thus depends on all key pair generation,
    signature generation, and signature verification files. That file
    also contains some simple performance benchmarks.

  - All other files support some utility code which is used by all three
    scheme components.

To keep the API simple, keys are referenced only in encoded format,
as byte slices. For a given degree (normally 512 or 1024), signing key,
verifying key, and signatures have a fixed size.

Key pair generation and signature verification use only plain integer code,
with no floating-point operations. Signature generation uses floating-point:

  - On some specific architectures, the native hardware is known to
    behave properly, i.e. exactly as specified as IEEE-754 for the
    "binary64" type with no extra precision or double-rounding issues. A
    few assembly instructions are furthermore used to ensure that some
    operations (rounding and square roots, mostly) are done with the raw
    hardware opcodes, avoiding any constant-time trouble with library
    functions. These optimized architectures are: 32-bit x86 with SSE2
    support (386.sse2), 64-bit x86 (amd64), 64-bit ARM (arm64) and
    64-bit RISC-V (riscv64). Note that when compiling for 32-bit x86,
    Go's defaults include SSE2 support, so that this case is optimized.

  - On all other architectures, integer-based emulation code is used.
    This code is presumed constant-time (to the extent that
    constant-time efforts are not ruined by some overly smart compiler
    or in-CPU JIT layer). The emulated code is about 8 times slower than
    the native code, hence not a record breaker, but not unusable
    either.

  - The integer-based emulation code can be selected at compile-time,
    even on architectures where hardware is supported, by defining the
    tag `fndsa_fp_emu`. This might be convenient for test purposes.
