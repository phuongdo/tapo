package org.cnrs.crbm.lib.io;

public class Row {

	String protein;
	int length;
	int byEyeTot;
    int byEyeTotSimply;
    int conclusion;
    int pdbBegin;
    int pdbEnd;
    int n;


    public int getByEyeTotSimply() {
        return byEyeTotSimply;
    }

    public void setByEyeTotSimply(int byEyeTotSimply) {
        this.byEyeTotSimply = byEyeTotSimply;
    }


    public int getPdbBegin() {
        return pdbBegin;
    }

    public void setPdbBegin(int pdbBegin) {
        this.pdbBegin = pdbBegin;
    }

    public int getPdbEnd() {
        return pdbEnd;
    }

    public void setPdbEnd(int pdbEnd) {
        this.pdbEnd = pdbEnd;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public int getConclusion() {
        return conclusion;
    }

	public void setConclusion(int conclusion) {
		this.conclusion = conclusion;
	}

	public String getProtein() {
		return protein;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public int getByEyeTot() {
		return byEyeTot;
	}

	public void setByEyeTot(int byEyeTot) {
		this.byEyeTot = byEyeTot;
	}

	public String getPdbCode() {

		return this.protein.substring(0, 4);
	}

	public String getPdbChain() {

		return this.protein.substring(5, 6);
	}

	@Override
	public String toString() {

        return protein + "\t" + length + "\t" + this.byEyeTot + "\t" + pdbBegin + "\t" + pdbEnd + "\t" + n + "\t" + byEyeTotSimply;
    }
}