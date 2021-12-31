package org.cnrs.crbm.lib.io;

import java.util.Date;

public class RowRepeatDB {

	String entry;
	String annlevel;
	String strclass;

	String region;
	String units;
	String insertions;

	public String getEntry() {
		return entry;
	}

	public void setEntry(String entry) {
		this.entry = entry;
	}

	public String getAnnlevel() {
		return annlevel;
	}

	public void setAnnlevel(String annlevel) {
		this.annlevel = annlevel;
	}

	public String getRegion() {
		return region;
	}

	public void setRegion(String region) {
		this.region = region;
	}

	public String getUnits() {
		return units;
	}

	public void setUnits(String units) {
		this.units = units;
	}

	public String getInsertions() {
		return insertions;
	}

	public void setInsertions(String insertions) {
		this.insertions = insertions;
	}

	public String getStrclass() {
		return strclass;
	}

	public void setStrclass(String strclass) {
		this.strclass = strclass;
	}

	public String getPdbCode() {

		return this.entry.substring(0, 4);
	}

	public String getPdbChain() {

		return this.entry.substring(4, 5);
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return entry + "\t" + annlevel;
	}

	@Override
	public boolean equals(Object obj) {

		RowRepeatDB row = (RowRepeatDB) obj;
		return super.equals(this.entry.equals(row.getEntry()));
	}
}