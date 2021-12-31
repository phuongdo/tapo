package org.cnrs.crbm.lib.classification;

import java.util.ArrayList;
import java.util.List;

public class ProData {

	String protein;
	String strClass;
	String anno;
	int nrTRs;
	int len;
	List<String> domains = new ArrayList<String>();

	public String getProtein() {
		return protein;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}

	public String getStrClass() {
		return strClass;
	}

	public void setStrClass(String strClass) {
		this.strClass = strClass;
	}

	public String getAnno() {
		return anno;
	}

	public void setAnno(String anno) {
		this.anno = anno;
	}

	public int getNrTRs() {
		return nrTRs;
	}

	public void setNrTRs(int nrTRs) {
		this.nrTRs = nrTRs;
	}

	public int getLen() {
		return len;
	}

	public void setLen(int len) {
		this.len = len;
	}

	public List<String> getDomains() {
		return domains;
	}

	public void setDomains(List<String> domains) {
		this.domains = domains;
	}

}
