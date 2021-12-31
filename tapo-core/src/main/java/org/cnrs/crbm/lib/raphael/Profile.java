package org.cnrs.crbm.lib.raphael;

public class Profile {

	double totalscore = 0.0;
	int[] peaks = null;
	int[] p = null;

	public double getTotalscore() {
		return totalscore;
	}

	public Profile(double totalscore, int[] p) {
		super();
		this.totalscore = totalscore;
		this.p = p;
	}

	public void setTotalscore(double totalscore) {
		this.totalscore = totalscore;
	}

	public int[] getP() {
		return p;
	}

	public void setP(int[] p) {
		this.p = p;
	}

}
