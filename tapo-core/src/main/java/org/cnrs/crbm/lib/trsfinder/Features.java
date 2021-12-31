package org.cnrs.crbm.lib.trsfinder;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.cnrs.crbm.lib.repeats.module.ProVector;

public class Features {

	String pdbCode;
	String pdbChain;
	String output;
	String strSS;
	String strAcc;
	String strSeqAp;
	String strSeq;
	String strSeqSADB;
	Atom[] atoms;
	List<ProVector> vectors = new ArrayList<ProVector>();
	double[] cmHistos;
	Chain chain;

	public void setCmapsList(Map<Integer, List<Integer>> cmapsList) {
		this.cmapsList = cmapsList;
	}

	public Map<Integer, List<Integer>> getCmapsList() {
		return cmapsList;
	}

	private Map<Integer, List<Integer>> cmapsList = new HashMap<Integer, List<Integer>>();


	public Features(String pdbCode, String pdbChain, String output,
			String strSS, String strAcc, String strSeqAp, String strSeq,
			Atom[] atoms, List<ProVector> vectors, double[] cmHistos) {
		super();
		this.pdbCode = pdbCode;
		this.pdbChain = pdbChain;
		this.output = output;
		this.strSS = strSS;
		this.strAcc = strAcc;
		this.strSeqAp = strSeqAp;
		this.strSeq = strSeq;
		this.atoms = atoms;
		this.vectors = vectors;
		this.cmHistos = cmHistos;
	}

	public Features() {

	}

	public String getPdbCode() {
		return pdbCode;
	}

	public void setPdbCode(String pdbCode) {
		this.pdbCode = pdbCode;
	}

	public String getPdbChain() {
		return pdbChain;
	}

	public void setPdbChain(String pdbChain) {
		this.pdbChain = pdbChain;
	}

	public String getOutput() {
		return output;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public String getStrSS() {
		return strSS;
	}

	public void setStrSS(String strSS) {
		this.strSS = strSS;
	}

	public String getStrAcc() {
		return strAcc;
	}

	public void setStrAcc(String strAcc) {
		this.strAcc = strAcc;
	}

	public String getStrSeqAp() {
		return strSeqAp;
	}

	public void setStrSeqAp(String strSeqAp) {
		this.strSeqAp = strSeqAp;
	}

	public String getStrSeq() {
		return strSeq;
	}

	public void setStrSeq(String strSeq) {
		this.strSeq = strSeq;
	}

	public Atom[] getAtoms() {
		return atoms;
	}

	public void setAtoms(Atom[] atoms) {
		this.atoms = atoms;
	}

	public List<ProVector> getVectors() {
		return vectors;
	}

	public void setVectors(List<ProVector> vectors) {
		this.vectors = vectors;
	}

	public double[] getCmHistos() {
		return cmHistos;
	}

	public void setCmHistos(double[] cmHistos) {
		this.cmHistos = cmHistos;
	}

	public String getStrSeqSADB() {
		return strSeqSADB;
	}

	public void setStrSeqSADB(String strSeqSADB) {
		this.strSeqSADB = strSeqSADB;
	}
	public Chain getChain() {
		return chain;
	}

	public void setChain(Chain chain) {
		this.chain = chain;
	}


}
