Command-line Fast Fourier Transform utility

This tool implements the standard FFT algorithm as a standalone
command-line utility.  The intention is to fit into the standard
linux/unix shell paradigm:  standard shell operations (piping inputs
and outputs, file etc.), character streams for input and output, etc.

The fft uses a mixed-radix algorithm, so the inputs are not restricted
to lengths that are a power of two.  The input length is broken into
its prime factors, and the fft uses the prime factors as the "fan out"
at each recursive step.  The fft algorithm is fastest if the largest
prime factor of N is small!  In the degenerate case that N is prime,
the algorithm runs in quadratic time.  If your dataset has a length
N that happens to have large prime factors, it is worth considering
trimming one or two elements off of the dataset.  A number near N
is likely to have a very different configuration of prime factors,
with a smaller maximum prime factor.
In that case, the resulting FFT is near the original, and it will be computed
much faster.

In some cases, use of this fft tool can avoid the need to use windows (Blackman, Hamming, etc.).
A standard use of windows is to minimize leakage around the FFT output frequencies due to
the fact that the sample length is not a multiple of the fundamental frequency of the data.

Take the FFT of the original dataset, and deduce the (approximate) fundamental
frequency of the data.  Then trim the data a little bit, so that the length is
an even multiple of the fundamental frequency.  The FFT of the modified dataset may have
less leakage around the calculated frequencies, because it is the correct length relative to
its fundamental frequency.

Here is a simple example:

    fft
    1
    0
    1
    0
    1
    0
    <eof>
     1.224744871392 0.000000000000
     0.000000000000 0.000000000000
    -0.000000000000 0.000000000000
     1.224744871392 0.000000000000
     0.000000000000 0.000000000000
    -0.000000000000 0.000000000000

The default input is double-precision real values, and
the default output is cartesian complex numbers (X + iY).
These input and output formats can be overridden with command-line arguments.

The following operation represents the identity function on
its inputs.  The second invocation generates the inverse FFT, and
takes cartesian complex numbers as input:

    fft | fft -b -ic -oa
    1
    2
    3
    101
    105
    <eof>
    1.000000000000
    2.000000000000
    2.999999999999
    101.000000000001
    105.000000000000

Here are the command-line arguments:

    usage:  fft
            [-h]         # help (this message)
            [-check N]   # calculate largest prime factor of N; smaller value => faster fft
            [-b]         # inverse (backward) fft

                         # default input format is real (non-complex) values
            [-ip]        # input complex values in polar (radius, radians)
            [-ic]        # input complex values in cartesian (real, imaginary)

                         # default output format is cartesian (real, imaginary) values
            [-op]        # output complex values in polar (radius, radians)
            [-oa]        # output double-precision amplitudes

            [input_file] # default input is stdin
