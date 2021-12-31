package org.cnrs.crbm.lib.math;

public class Threshold {

	// f(x) = MaxRMSD/(1+e(-x))
	/**
	 * 
	 * @param len
	 *            of alignment
	 * @param maxRMSD
	 * @return
	 */
	static public double logistic(double len, double maxRMSD) {

		return (maxRMSD / (1 + (Math.exp(-len))));

	}

	/**
	 * FIXED = 4 Ao
	 * 
	 * @param x

     * @return
     */
    static public double functRMSD(double x) {

		// return 4;
        //return 6.0;
        return funct2(x, 5.2);
        // return x / Math.pow(x, 0.6);
    }

	static public double functGapsRes(double x) {
		if (x < 18)
			return 0;
        else {

            return 1.889e-01 * x - 1.256e-15;

        }
    }

	/**
	 * aL = Ra/R
	 * 
	 * @param x
	 *            number of residues in alignment
	 * @return aL
	 */
	static public double functaL(double x) {

		return x * 0.71;

	}

	/**
	 * aL = Ra/S
	 * 
	 * @param x
	 *            number of residues in alignment
	 * @return aL
	 */
	static public double functaS(double x) {

		return x * 0.86;

	}

	static public double funct2(double x, double max) {
		if (x < 25)
			return (max * x) / (max + 3 * x);
		else if (x < 30)
			return (max * x + x) / (max + x + x);
		else
			return max * x * 2 / (max + x + x);
	}

	/**
	 * y = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5
	 * 
	 * Coefficients: p1 = 0.0012812 p2 = -0.12044 p3 = 3.9703 p4 = -52.258 p5 =
	 * 264.6
	 * 
	 * 
	 * 
	 * @param x
	 * @return
	 */
	static public double funcThresholdForSA(double x) {

		return funCubic(x);
	}

	/**
	 * y = p1*x^2 + p2*x + p3
	 * 
	 * Coefficients: p1 = -0.025698 p2 = 2.3949 p3 = 3.2385
	 * 
	 * Norm of residuals = 39.827
	 * 
	 * @param x
	 * @return
	 */
	static public double funcQuadraticSA(double x) {
		double p1 = -0.025698;
		double p2 = 2.3949;
		double p3 = 3.2385;
		return p1 * Math.pow(x, 2) + p2 * x + p3;

	}

	/**
	 * y = p1*x^3 + p2*x^2 + p3*x + p4
	 * 
	 * Coefficients: p1 = -0.0011053 p2 = 0.1927 p3 = -7.1171 p4 = 87.518
	 * 
	 * Norm of residuals = 27.065
	 * 
	 * @param x
	 * @return
	 */
	static public double funCubic(double x) {
		double p1 = -0.0011053;
		double p2 = 0.1927;
		double p3 = -7.1171;
		double p4 = 87.518;
		return p1 * Math.pow(x, 3) + p2 * Math.pow(x, 2) + p3 * x + p4;

	}


    static public double funcTmScore(double x) {
        double p1 = 2.289e-08;
        double p2 = -4.840e-05;
        double p3 = 1.329e-02;
        double p4 = -1.133e-01;
        double y = p1 * Math.pow(x, 3) + p2 * Math.pow(x, 2) + p3 * x + p4;

        if (x > 100) {
            return 0.8;
        } else {
            if (y < 0.17) return 0.17;
            else
                return y;

        }

    }

	public static void main(String[] args) {
		for (int i = 0; i < 200; i++) {
            System.out.println(i + "\t" + Threshold.funcTmScore(i));
        }

	}
}
