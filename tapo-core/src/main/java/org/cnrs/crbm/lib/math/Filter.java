package org.cnrs.crbm.lib.math;

import org.jtransforms.fft.DoubleFFT_1D;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Filter {

    public static double[] filter(double[] v) throws Exception {
        // using gausse filter

        //double[] f = new double[v.length];

        int n = v.length;
        int mu = 6; // 6: default value is 20
        int m = 40; // 40: total size of the filter, odd number is better default value is 101

        double[] h = buildGaussianFilter(m, (double) mu / (4 * n), n);
        return convolve(v, h);

    }


    private static double sqr(double x) {
        return x * x;
    }

    public static double[] fftAutoCorrelation(double[] magFFT) {
        int magCnt = magFFT.length;
        DoubleFFT_1D fft = new DoubleFFT_1D(magCnt);
        fft.realForward(magFFT);

        magFFT[0] = (magFFT[0] * magFFT[0]);
        for (int i = 1; i < (magCnt - (magCnt % 2)) / 2; i++) {
            magFFT[2 * i] = magFFT[2 * i] * magFFT[2 * i] + magFFT[2 * i + 1] * magFFT[2 * i + 1];
            magFFT[2 * i + 1] = 0.0;
        }

        if (magCnt % 2 == 0) {
            magFFT[1] = (magFFT[1] * magFFT[1]);
        } else {
            magFFT[magCnt / 2] = (magFFT[magCnt - 1] * magFFT[magCnt - 1] + magFFT[1] * magFFT[1]);
        }

        double[] autocorr = new double[magCnt];
        System.arraycopy(magFFT, 0, autocorr, 0, magCnt);
        DoubleFFT_1D ifft = new DoubleFFT_1D(magCnt);
        ifft.realInverse(autocorr, false);

        for (int i = 1; i < autocorr.length; i++)
            autocorr[i] /= autocorr[0];
        autocorr[0] = 1.0;

        return autocorr;
    }

    public static double[] fastFilter(double[] v) throws Exception {
        double[] f = fastFilter3(v);
        f = fastFilter6(f);
        return f;
    }

    public static double[] fastFilter3(double[] v) throws Exception {
        double[] f = new double[v.length];
        f[0] = v[0];
        for (int i = 1; i < v.length - 1; i++) {
            f[i] = (double) (v[i - 1] + v[i] + v[i + 1]) / 3;
        }
        f[v.length - 1] = v[v.length - 1];
        return f;
    }

    public static double[] fastFilter6(double[] v) throws Exception {
        double[] f = new double[v.length];
        int winsize = 6;
        for (int i = 1; i < v.length - winsize; i++) {

            for (int j = 0; j < winsize; j++) {
                f[i] += (double) v[i + j];
            }

            f[i] = (double) f[i] / winsize;
        }

        for (int i = v.length - winsize; i < v.length; i++) {
            f[i] = f[v.length - winsize - 1];
        }

        return f;
    }

    public static double[] buildGaussianFilter(int m, double sigma, int n) {

        // n = m;
        //
        // x = ( (0:n-1)-(n-1)/2 )/(N-1);
        // f = exp( -x.^2/(2*s^2) );
        // f = f / sum(f(:));
        double[] x = new double[n];
        for (int i = 0; i < n; i++)
            x[i] = (double) (i - (n - 1) / 2) / (n - 1);

        double[] f = new double[n];
        for (int i = 0; i < n; i++) {
            f[i] = Math.exp(-Math.pow(x[i], 2) / (2 * Math.pow(sigma, 2)));

        }

        double sum = sum(f);

        for (int i = 0; i < n; i++)
            f[i] = (double) f[i] / sum;
        return f;

    }

    public static double[] convolve(double[] v, double[] h) throws Exception {

        int n = v.length;
        int p = h.length;

        if (p > n)
            throw new Exception("h filter should be shorter than x.");

        double[] y = new double[n];

        int d1 = (int) Math.floor((p) / 2);// padding before
        int d2 = p - d1 - 1; // padding after

        // double[] xx = new double[d1 + v.length + d2];
        // xx = [ x(d1:-1:1); x; x(end:-1:end-d2+1) ];
        List<Double> xx = new ArrayList<Double>();

        for (int i = d1 - 1; i > 0; i--) {

            xx.add(v[i]);
        }
        for (int i = 0; i < v.length; i++) {

            xx.add(v[i]);
        }
        for (int i = v.length - 1; i > v.length - d2; i--) {
            xx.add(v[i]);
        }

        // y = conv(xx,h);

        y = conv(xx.toArray(new Double[xx.size()]), h);

        // y = y(p:end-p+1);

        List<Double> r = new ArrayList<Double>();

        for (int i = p; i < y.length - p + 4; i++) {
            r.add(y[i]);

        }

        double[] k = new double[r.size()];

        int i = 0;
        for (Double t : r) {

            k[i] = t;
            i++;
        }
        return k;

        //
        // h = [ h(d+1:end); zeros(n-p,1); h(1:d) ];
        // y = real( ifft( fft(x).*fft(h) ) );

    }

    public static double[] conv(Double[] v, double[] h) {

        int SignalLen = v.length;
        int KernelLen = h.length;
        double[] c = new double[SignalLen + KernelLen - 1];

        for (int n = 0; n < SignalLen + KernelLen - 1; n++) {
            int kmin, kmax, k;

            c[n] = 0;

            kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
            kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

            for (k = kmin; k <= kmax; k++) {
                c[n] += v[k] * h[n - k];
            }
        }

        return c;

    }

    public static double sum(double[] v) {

        double sum = 0.0;
        for (int i = 0; i < v.length; i++)
            sum += v[i];
        return sum;
    }

}
