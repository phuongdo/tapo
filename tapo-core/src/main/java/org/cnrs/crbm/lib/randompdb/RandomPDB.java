package org.cnrs.crbm.lib.randompdb;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.random.RandomDataGenerator;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.jama.Matrix;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

public class RandomPDB {

	List<Atom> feedAtoms = new ArrayList<Atom>();

	public static void main(String[] args) throws Exception {

		RandomPDB random = new RandomPDB();

		random.build();

	}

	public String build() throws Exception {
		Random random = new Random();
		List<String> pdbs = new ArrayList<String>();
		pdbs.add("2KIV");
		pdbs.add("3A1Y");
		pdbs.add("4IRT");
		pdbs.add("1A4Y");
		pdbs.add("1LXA");
		pdbs.add("1R56");
		pdbs.add("1A17");
		pdbs.add("1IA5");
		pdbs.add("1GVM");
		List<List<Atom>> listSS = new ArrayList<List<Atom>>();
		// get template
		for (String pdbCode : pdbs) {
			listSS.addAll(this.generate(pdbCode));
		}

		// build a protein;

		int noElement = 10;

		List<Atom> fakeProtein = new ArrayList<Atom>();
		for (int i = 0; i < noElement; i++) {
			// get random from listSS
			try {
				List<Atom> ss = listSS.get(random.nextInt(listSS.size()));
				// for (int t = 0; t < 10; t++)
				// rotateAtoms(ss);
				// build loop;
				List<Atom> loop = this.randomWalk(random.nextInt(20) + 10);
				// rotateAtoms(loop);
				if (fakeProtein.size() == 0)
					fakeProtein = join(ss, loop);
				else
					fakeProtein = join(fakeProtein, join(ss, loop));
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

		// String AAs = "AVILMFYWSTNQCUGPRHKDE";
		StringBuilder builder = new StringBuilder();
		int i = 1;
		for (Atom atom : fakeProtein) {

			int resSeq = i;

			Atom feed = feedAtoms.get(random.nextInt(feedAtoms.size() - 1));
			String resName = feed.getGroup().getPDBName();
			feed.setCoords(atom.getCoords());
			builder.append(this.getPDBLine(feed, i, resName, resSeq) + "\n");
			i++;
		}

		// write to pdb file
		// String outputfile = "output/fake.pdb";
		// FileOutputStream out = new FileOutputStream(outputfile);
		// PrintStream p = new PrintStream(out);
		// p.write(builder.toString().getBytes());

		// save show

		// this.saveShowCommand(fakeProtein);
		return builder.toString();

	}

	private void saveShowCommand(List<Atom> atoms) {

		String outfile_pml = "output/fake_show.pml";

		try {
			FileWriter fstream_pml = new FileWriter(outfile_pml);
			BufferedWriter out_pml = new BufferedWriter(fstream_pml);

			// Structure newstruc = new StructureImpl();
			// Chain c1 = new ChainImpl();
			// c1.setChainID(pdbChain);
			out_pml.write("set dash_gap, 0.0\n");

			int j = 0;
			for (int i = 1; i < atoms.size(); i++) {

				// System.out.println(v.getStart() + ":"
				// + v.getOriginalStart());
				int p1 = i;
				int p2 = i + 1;
				out_pml.write("distance d" + i + ", ////" + p1 + ", ////" + p2
						+ " \n");
				out_pml.write("hide labels, d" + i + "\n");
				j++;

			}

			out_pml.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private List<Atom> join(List<Atom> ss, List<Atom> loop) {
		List<Atom> join = new ArrayList<Atom>();

		Atom vector_shift = new AtomImpl();
		vector_shift.setX(ss.get(ss.size() - 1).getX() - loop.get(0).getX());
		vector_shift.setY(ss.get(ss.size() - 1).getY() - loop.get(0).getY());
		vector_shift.setZ(ss.get(ss.size() - 1).getZ() - loop.get(0).getZ());

		// move loop to near ss

		// for (Atom atom : loop) {
		// Calc.shift(atom, vector_shift);
		// }

		join.addAll(ss);
		for (int i = 1; i < loop.size(); i++) {
			Calc.shift(loop.get(i), vector_shift);
			join.add(loop.get(i));
		}

		return join;
	}

	private List<List<Atom>> generate(String pdbCode) {
		List<List<Atom>> listSS = new ArrayList<List<Atom>>();

		try {
			String pdbChain = "A";
			// get data
			// String pdbFile = PdbTools.downloadPDB(pdbCode);
			Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
			Atom[] atoms = StructureTools.getAtomCAArray(structure
					.getChainByPDB(pdbChain));

			for (int i = 0; i < atoms.length; i++)
				feedAtoms.add(atoms[i]);
			DSSP dssp = new DSSP(pdbCode);
			// make a filter first
			String secondaryStr = dssp.filterSS(dssp.getCombinedSS(atoms,
					pdbChain));
			List<ProVector> vectors = VectorShape.getVectors(atoms,
					secondaryStr);

			vectors = ProVector.toSecondaryVector(vectors);
			// 1 second
			// for (ProVector vector : vectors) {
			//
			// List<Atom> latom = new ArrayList<Atom>();
			// for (int i = vector.getPosStart(); i <= vector.getPosEnd(); i++)
			// {
			// latom.add((Atom) atoms[i].clone());
			// }
			//
			// listSS.add(latom);
			// }

			int winsize = 2;
			for (int i = 0; i < vectors.size() - winsize + 1; i++) {
				List<ProVector> sub_list = new ArrayList<ProVector>();
				for (int j = 0; j < winsize; j++) {
					ProVector v = vectors.get(i + j);
					sub_list.add(v);
				}
				List<Atom> latom = new ArrayList<Atom>();
				for (int k = sub_list.get(0).getPosStart(); k <= sub_list.get(
						sub_list.size() - 1).getPosEnd(); k++) {
					latom.add((Atom) atoms[k].clone());
				}
				listSS.add(latom);

			}
			winsize = 3;
			for (int i = 0; i < vectors.size() - winsize + 1; i++) {
				List<ProVector> sub_list = new ArrayList<ProVector>();
				for (int j = 0; j < winsize; j++) {
					ProVector v = vectors.get(i + j);
					sub_list.add(v);
				}
				List<Atom> latom = new ArrayList<Atom>();
				for (int k = sub_list.get(0).getPosStart(); k <= sub_list.get(
						sub_list.size() - 1).getPosEnd(); k++) {
					latom.add((Atom) atoms[k].clone());
				}
				listSS.add(latom);

			}

			// randomWalk(atoms);
		} catch (StructureException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return listSS;

	}

	private List<Atom> randomWalk(int NC) throws Exception {

		Random random = new Random();
		List<Atom> loop = new ArrayList<Atom>();

		// int NC = 10;
		RandomDataGenerator randomData = new RandomDataGenerator();
		double NN = 5.5, TC = 4.0, F = 0.5;
		double lc = TC - F;
		double uc = TC + F;
		double ln = NN - F;
		double un = NN + F;
		double[][] C = new double[NC][3];
		C[0][0] = 0.0;
		C[0][1] = 0.0;
		C[0][2] = 0.0;
		C[1][0] = 4.0;
		C[1][1] = 0.0;
		C[1][2] = 0.0;

		for (int i = 2; i < NC; i++) {
			int j = 1, t = 0, n = -3;
			double[] X = new double[3];
			do {
				// Loop to try for an accepted point
				t++;

				X[0] = C[i - 1][0] - uc + (2 * uc * Math.random());
				X[1] = C[i - 1][1] - uc + (2 * uc * Math.random());
				X[2] = C[i - 1][2] - uc + (2 * uc * Math.random());
				double v = dist(X, C[i - 1]);

				if (v < lc || v > uc)
					continue;

				for (int k = 0; k < i - 1; k++) {
					v = dist(X, C[k]);
					if (v < ln) {
						j = 1;
						break;
					}
					if (v < un) {
						j--;
					}
				}

				if (j < n) {
					j = 0;
				} else {
					j = 1;
					if (t > 1000) {
						t = 0;
						n++;
					}
				}

			} while (j == 1 && i > 2 && n <= 0);

			C[i][0] = X[0];
			C[i][1] = X[1];
			C[i][2] = X[2];
		}

		for (int i = 0; i < NC; i++) {

			Atom atom = new AtomImpl();
			atom.setCoords(C[i]);

			loop.add(atom);

		}

		return loop;
		// translate coordinate to atoms
		// Atom[] randomAtoms = new Atom[NC];
		// create new Structure
		// Structure struc = new StructureImpl();
		// Chain c = new ChainImpl();
		// StringBuilder builder = new StringBuilder();
		//
		// for (int i = 0; i < NC; i++) {
		// // get atom from random;
		// int rnumber = random.nextInt(atoms.length);
		// // Atom atom = (Atom) atoms[random.nextInt(atoms.length)].clone();
		// Atom atom = (Atom) atoms[rnumber].clone();
		// // Atom atom = atoms[i];
		// atom.setGroup(atoms[rnumber].getGroup());
		//
		// // double[] vector_tran = new double[3];
		// Atom vector_shift = new AtomImpl();
		// vector_shift.setX(C[i][0] - atom.getX());
		// vector_shift.setY(C[i][1] - atom.getY());
		// vector_shift.setZ(C[i][2] - atom.getZ());
		// // atom.setCoords(C[i]);
		//
		// AminoAcid aa = (AminoAcid) atom.getGroup();
		// String resName = aa.getPDBName();
		// int resSeq = i + 1;
		//
		// // g.getna
		// int j = i + 1;
		//
		// Calc.shift(atom.getGroup(), vector_shift);
		// for (Atom a : aa.getAtoms()) {
		//
		// builder.append(this.getPDBLine(a, i, resName, resSeq) + "\n");
		// j++;
		// }
		//
		// // System.out.println(atom);
		// // atom.getGroup().clearAtoms();
		// // atom.getGroup().addAtom(atom);
		// // ResidueNumber resNum = new ResidueNumber();
		// // resNum.setSeqNum(i);
		//
		// // atom.getGroup().setResidueNumber(resNum);
		//
		// // c.addGroup(atom.getGroup());
		// // randomAtoms[i] = atom;
		//
		// }

		// write to pdb file
		// String outputfile = "output/fake.pdb";
		// FileOutputStream out = new FileOutputStream(outputfile);
		// PrintStream p = new PrintStream(out);
		// p.write(builder.toString().getBytes());

	}

	public void rotateAtoms(List<Atom> ca) {
		Matrix rotMatrix = this.genRandomRot();
		for (int i = 0; i < ca.size(); i++) {
			Calc.rotate(ca.get(i), rotMatrix);
		}
	}

	public Matrix genRandomRot() {
		Matrix rotMatrix = new Matrix(3, 3);
		double rangeMin = 0;
		double rangeMax = 360;
		Random random = new Random();

		double ang_x = rangeMin + (rangeMax - rangeMin) * random.nextDouble();
		double ang_y = rangeMin + (rangeMax - rangeMin) * random.nextDouble();
		double ang_z = rangeMin + (rangeMax - rangeMin) * random.nextDouble();
		Matrix r_x = new Matrix(3, 3);
		r_x.set(1, 1, Math.cos(ang_x));
		r_x.set(1, 2, -Math.sin(ang_x));
		r_x.set(2, 2, Math.cos(ang_x));
		r_x.set(1, 2, Math.sin(ang_x));
		r_x.set(0, 0, 1);
		Matrix r_y = new Matrix(3, 3);
		r_y.set(0, 0, Math.cos(ang_y));
		r_y.set(0, 2, -Math.sin(ang_y));
		r_y.set(2, 2, Math.cos(ang_y));
		r_y.set(2, 0, Math.sin(ang_y));
		r_y.set(1, 1, 1);

		Matrix r_z = new Matrix(3, 3);
		r_z.set(0, 0, Math.cos(ang_z));
		r_z.set(0, 1, -Math.sin(ang_z));
		r_z.set(1, 1, Math.cos(ang_z));
		r_z.set(1, 0, Math.sin(ang_z));
		r_z.set(2, 2, 1);

		rotMatrix = r_x.times(r_y.times(r_z));

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				rotMatrix.set(i, j, random.nextDouble() - 1);
			}
		}

		return rotMatrix;
	}

	private String getPDBLine(Atom a, int oderOfAtom, String resName, int resSeq) {
		String formatted = String.format(
				"ATOM%7s  %-4s%-3s %-1s%4s%12s%8s%8s%6s%6s           %-3s",
				oderOfAtom, a.getName(), resName, "A", resSeq,
				NumberFormatUtils.format(a.getX()),
				NumberFormatUtils.format(a.getY()),
				NumberFormatUtils.format(a.getZ()), a.getOccupancy(),
				a.getTempFactor(), a.getElement());

		return formatted;
	}

	private double dist(double[] X, double[] C) {

		return Math.sqrt((X[0] - C[0]) * (X[0] - C[0]) + (X[1] - C[1])
				* (X[1] - C[1]) + (X[2] - C[2]) * (X[2] - C[2]));

	}
}
