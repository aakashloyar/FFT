package com.aakash.Polynomail_Multiplication;

public class FFT
{
    public static class Complex
    {
        public final double re; // real part
        public final double im; // imaginary part


        public Complex(double real, double imag) {
            re = real;
            im = imag;
        }

        public String toString() {
            if (im == 0)
                return re + "";
            if (re == 0)
                return im + "i";
            if (im < 0)
                return re + " - " + (-im) + "i";
            return re + " + " + im + "i";
        }
        public Complex times(double alpha) {
            return new Complex(alpha * re, alpha * im);
        }

        public Complex conjugate() {
            return new Complex(re, -im);
        }

        public double abs() {
            return Math.sqrt((Math.pow(re,2)+Math.pow(im,2)));
        }

        public double phase() {
            return Math.atan2(im, re);
        }

        public Complex plus(Complex b) {
            Complex a = this;
            double real = a.re + b.re;
            double imag = a.im + b.im;
            return new Complex(real, imag);
        }

        public Complex minus(Complex b) {
            Complex a = this;
            double real = a.re - b.re;
            double imag = a.im - b.im;
            return new Complex(real, imag);
        }

        public Complex times(Complex b) {
            Complex a = this;
            double real = a.re * b.re - a.im * b.im;
            double imag = a.re * b.im + a.im * b.re;
            return new Complex(real, imag);
        }

        public Complex reciprocal() {
            double scale = re * re + im * im;
            return new Complex(re / scale, -im / scale);
        }

        public Complex sin() {
            return new Complex(Math.sin(re) * Math.cosh(im), Math.cos(re)
                    * Math.sinh(im));
        }


        public Complex cos() {
            return new Complex(Math.cos(re) * Math.cosh(im), -Math.sin(re)
                    * Math.sinh(im));
        }


        public Complex tan() {
            return sin().divides(cos());
        }
        public Complex divides(Complex b) {
            Complex a = this;
            return a.times(b.reciprocal());
        }


        public Complex exp() {
            return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re)
                    * Math.sin(im));
        }

        public static Complex[] fft(Complex[] x) {
            int N = x.length;


            if (N == 1)
                return new Complex[] { x[0] };


            if (N % 2 != 0) throw new RuntimeException("N is not a power of 2");
            // fft of even terms
            Complex[] even = new Complex[N / 2];
            for (int k = 0; k < N / 2; k++) {
                even[k] = x[2 * k];
            }
            Complex[] q = fft(even);

            // fft of odd terms
            Complex[] odd = even;
            for (int k = 0; k < N / 2; k++) {
                odd[k] = x[2 * k + 1];
            }
            Complex[] r = fft(odd);


            Complex[] y = new Complex[N];
            for (int k = 0; k < N / 2; k++) {
                double kth = -2 * k * Math.PI / N;
                Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
                y[k] = q[k].plus(wk.times(r[k]));
                y[k + N / 2] = q[k].minus(wk.times(r[k]));
            }
            return y;
        }


        public static Complex[] ifft(Complex[] x) {
            int N = x.length;
            Complex[] y = new Complex[N];


            for (int i = 0; i < N; i++) {
                y[i] = x[i].conjugate();
            }


            y = fft(y);


            for (int i = 0; i < N; i++) {
                y[i] = y[i].conjugate();
            }


            for (int i = 0; i < N; i++) {
                y[i] = y[i].times(1.0 / N);
            }

            return y;

        }


        public static void print(Complex[] arr, String s) {
            System.out.println(s);
            for (int i = 0; i < arr.length; i++) System.out.println(arr[i]);
            System.out.println();
        }

        public static void main(String[] args) {
            int N = 8;
            Complex[] p1 = new Complex[N];
            Complex[] p2 = new Complex[N];
            int[] arr1=new int[]{1,1,0,1,0,0,0,0};
            int[] arr2=new int[]{1,1,0,1,0,0,0,0};

            // original data
            System.out.println("Given Polynomial");
            for (int i = 0; i < N; i++) p1[i] = new Complex(arr1[i], 0);
            for (int i = 0; i < N; i++) p2[i] = new Complex(arr2[i], 0);
            print(p1, "p1");
            print(p2, "p2");


            //******STEP:1 : Polynomial to point form
            System.out.println("Fourier Transform of Given Polynomial");
            Complex[] y1 = fft(p1);
            print(y1, "y1 = fft(p1(x)");
            Complex[] y2 = fft(p2);
            print(y2, "y2 = fft(p2(x)");




            //******STEP:2 : Point form Multiplication
            System.out.println("Point multiplication form of given polynomial");
            Complex[] y = new Complex[N];
            for (int i = 0; i < N; i++)  y[i] = y1[i].times(y2[i]);

            print(y, "multiplied point form");



            //******STEP:3 : Taking the Inverse
            System.out.println("Inverse Fourier Transform of multiplied point form");
            Complex[] p = ifft(y);
            print(p, "p = ifft(y)");
            System.out.println("This is multiplied ploynomial");

        }

    }
}

