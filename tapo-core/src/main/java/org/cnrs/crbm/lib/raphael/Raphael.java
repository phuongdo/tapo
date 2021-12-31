package org.cnrs.crbm.lib.raphael;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.jama.Matrix;
import org.cnrs.crbm.lib.math.Filter;
import org.cnrs.crbm.lib.math.Peaks;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.*;

/**
 * <pre>
 * Methods paper: Ian Walsh, Francesco G. Sirocco, Giovanni Minervini, Carlo
 * Ferrari and Silvio C.E. Tosatto, Raphael: Recognition, periodicity and
 * residue-based assignment of structurally repeating proteins , IN PRESS.
 * (2012)
 * </pre>
 *
 * @author phuongdo
 */
public class Raphael {

    public static void main(String[] args) throws Exception {
        //String pdbCode = "1ia5";
//        String pdbCode = "1ozn";
        String pdbCode = "1mr7";
        String pdbChain = "A";
        // String pdbFile = PdbTools.downloadPDB(pdbCode);
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);

        Atom[] atoms = StructureTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain));

        System.out.println(atoms.length / 2);
        Atom[] atomsR = Fragement.getFragementsofAtoms(atoms, atoms.length / 2, atoms.length - 1);


        Raphael raphael = new Raphael(atomsR);
        System.out.println(raphael.getRepeatLength());
        System.out.println(raphael.getTotalScore());

    }

    private double repeatLength = 0.0;
    private double totalScore = 0.0;
    private double variancetScore = 0.0;
    private Atom[] atoms = null;
    private int T = 5;
    private int peakDetectorWindow = 5;
    private double stringency = 1;

    public Raphael() {
    }


    public Raphael(Atom[] atoms) {
        // clone
//        this.atoms = new Atom[atoms.length];
//        for (int i = 0; i < this.atoms.length; i++) {
//            this.atoms[i] = (Atom) atoms[i].clone();
//        }
        this.atoms = StructureTools.cloneAtomArray(atoms);
        try {
            run();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

//	public Raphael(Atom[] atoms, int T ) {
//		// clone
//		this.atoms = new Atom[atoms.length];
//		for (int i = 0; i < this.atoms.length; i++) {
//			this.atoms[i] = (Atom) atoms[i].clone();
//		}
//		this.T = T;
//
//		try {
//			run();
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	}

    public Raphael(List<Atom> atoms, int T, int peakDetectorWindow, double stringency) {
        // clone
        this.atoms = new Atom[atoms.size()];
        for (int i = 0; i < this.atoms.length; i++) {
            this.atoms[i] = (Atom) atoms.get(i).clone();
        }
        this.T = T;
        this.peakDetectorWindow = peakDetectorWindow;
        this.stringency = stringency;
        try {
            run();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private Atom[] cloneAtoms(Atom[] atoms) {
        // clone
        Atom[] clones = new Atom[atoms.length];
        for (int i = 0; i < this.atoms.length; i++) {
            clones[i] = (Atom) atoms[i].clone();
        }

        return clones;

    }

    private void run() throws Exception {
        DescriptiveStatistics stats_leng_finals = new DescriptiveStatistics();
        DescriptiveStatistics stats_score_finals = new DescriptiveStatistics();
        for (int step = 0; step < 30; step++) {
            List<Profile> plist = new ArrayList<Profile>();
            for (int i = 0; i < 200; i++) {
                Matrix rotMatrix = this.genRandomRot();
                Atom tranMatrix = this.genRandomTran();
                Atom[] clones = this.cloneAtoms(atoms);
                this.rotateShiftAtoms(clones, rotMatrix, tranMatrix);
                double[] px = getProfileX(clones);
                double[] py = getProfileY(clones);
                double[] pz = getProfileZ(clones);

                px = Filter.fastFilter(px);
                py = Filter.fastFilter(py);
                pz = Filter.fastFilter(pz);
//				this.show(px);
//				System.exit(0);
                // px = Filter.filter(px);
                // py = Filter.filter(pz);
                // pz = Filter.filter(py);

                Profile x = this.getProfile(px);
                Profile y = this.getProfile(py);
                Profile z = this.getProfile(pz);

                if (x.getTotalscore() > 0)
                    plist.add(x);
                if (y.getTotalscore() > 0)
                    plist.add(y);
                if (z.getTotalscore() > 0)
                    plist.add(z);

            }

            // get
            DescriptiveStatistics stats = new DescriptiveStatistics();
            for (int i = 0; i < atoms.length; i++) {
                // DescriptiveStatistics stats_i = new DescriptiveStatistics();
                int[] p = new int[plist.size()];
                for (int j = 0; j < plist.size(); j++) {
                    p[j] = plist.get(j).getP()[i];
                    // remove
                    // System.out.println(p[j]);
                    // stats_i.addValue(p[j]);

                    if (0 < p[j] && p[j] < 60) {
                        stats.addValue(p[j]);
                        // stats_i.addValue(p[j]);
                    }
                }
                // double mean_i = stats_i.getMean();
                // System.out.println(mean_i);
                // if (0 < mean_i && mean_i < 60)
                // stats.addValue(mean_i);
                // Map<Integer, Integer> F = getCount(p);
            }
            double P_avg = stats.getMean();
            double P_sd = stats.getStandardDeviation();

            stats = new DescriptiveStatistics();
            for (int i = 0; i < atoms.length; i++) {

                // DescriptiveStatistics stats_i = new DescriptiveStatistics();
                int[] p = new int[plist.size()];
                for (int j = 0; j < plist.size(); j++) {
                    p[j] = plist.get(j).getP()[i];
                    if (P_avg - P_sd / 2 <= p[j] && p[j] <= P_avg + P_sd / 2)
                        stats.addValue(p[j]);
                }
                // double mean_i = stats_i.getMean();
                // if (P_avg - P_sd / 2 <= mean_i && mean_i <= P_avg + P_sd / 2)
                // stats.addValue(mean_i);
                // Map<Integer, Integer> F = getCount(p);
            }

            DescriptiveStatistics stats_totalScore = new DescriptiveStatistics();
            for (int j = 0; j < plist.size(); j++) {
                stats_totalScore.addValue(plist.get(j).getTotalscore());
            }
            // System.out.println("Length: "
            // + NumberFormatUtils.format(stats.getMean()));

            if (stats_totalScore.getN() > 0) {
                stats_leng_finals.addValue(stats.getMean());
                stats_score_finals.addValue(stats_totalScore.getMean());
            }

        }

        if (stats_score_finals.getN() > 0) {
            this.repeatLength = stats_leng_finals.getMean();
            this.variancetScore = stats_leng_finals.getVariance();
            this.totalScore = stats_score_finals.getMean();

        }
        // System.out.println("sd:" + stats_leng_finals.getStandardDeviation());
    }

    public double getRepeatLength() {
        return repeatLength;

    }

    public double getTotalScore() {
        return this.totalScore;

    }

    public void show(double[] v) {

        for (int i = 0; i < v.length; i++)
            System.out.print(NumberFormatUtils.format(v[i]) + " ");
    }

    public void show(int[] v) {

        for (int i = 0; i < v.length; i++)
            System.out.print(NumberFormatUtils.format(v[i]) + " ");
    }

    public Profile getProfile(double[] inputStart) {


        //Peaks.findPeaks()
        //PeakDetector detector = new PeakDetector(inputStart);
        // get peaks
        //int[] peaks = detector.process(this.peakDetectorWindow, this.stringency);
        // get periods

        LinkedList<Integer> lstPeaks = Peaks.findPeaks(inputStart, 15, 0);
        int[] peaks = new int[lstPeaks.size()];
        for (int i = 0; i < peaks.length; i++) {
            peaks[i] = lstPeaks.get(i);
        }

        int[] profiles = new int[atoms.length];
        if (peaks.length > 1) {
            int[] periods = new int[peaks.length - 1];
            for (int i = 1; i < peaks.length - 1; i++) {

                periods[i] = peaks[i + 1] - peaks[i];
                // profiles[peaks[i]] = periods[i];
                for (int j = peaks[i]; j < peaks[i + 1]; j++) {
                    profiles[j] = periods[i];
                }
            }

            // for (int i = 0; i < peaks[0]; i++) {
            // profiles[i] = periods[0];
            // }
            // for (int i = peaks[peaks.length - 1]; i < atoms.length; i++) {
            // profiles[i] = periods[periods.length - 1];
            // }

            // get labels

            int[] labels = new int[periods.length];
            for (int i = 0; i < labels.length; i++)
                labels[i] = -1;

            int label = 0;
            for (int i = 0; i < periods.length; i++) {

                if (labels[i] == -1) {
                    int periodRf = periods[i];
                    labels[i] = label;

                    for (int j = i + 1; j < periods.length; j++) {
                        if (labels[j] == -1)
                            if (periodRf - this.T <= periods[j]
                                    && periods[j] <= periodRf + this.T) {
                                labels[j] = label;
                            }
                    }

                    label++;
                }

            }

            Map<Integer, Integer> C = this.getCount(labels);
            double p = 0.49;
            double W = this.getWindowScore(C, labels);
            double B = this.getBridgeScore(C, labels);
            int N = labels.length;// sequence length
            if (N > 1) {
                double totalscore = (p * W + (1 - p) * B) / ((N - 1) * N);
                return new Profile(totalscore, profiles);
            } else
                return new Profile(0, profiles);
        } else
            return new Profile(0, profiles);

    }

    public double getWindowScore(Map<Integer, Integer> C, int[] labels) {

        double score = 0.0;
        for (int i = 0; i < labels.length - 1; i++) {
            int l1 = labels[i];
            int l2 = labels[i + 1];
            if (l2 == l1)
                score += 2 * C.get(l1);
        }
        return score;
    }

    public double getBridgeScore(Map<Integer, Integer> C, int[] labels) {
        double score = 0.0;
        for (int i = 0; i < labels.length; i++) {
            int l1 = labels[i];
            int index = -1;
            for (int j = i + 1; j < labels.length; j++) {
                if (l1 == labels[j]) {
                    index = j;
                    break;
                }
            }
            if (index != -1 && (index - i) > 1) {// found
                double insertScore = 0;
                for (int k = i + 1; k < index; k++) {
                    insertScore += C.get(labels[k]);
                }
                score += 2 * C.get(l1) - insertScore;
            }
        }
        return score;
    }

    public Map<Integer, Integer> getCount(int[] labels) {
        Map<Integer, Integer> count = new HashMap<Integer, Integer>();
        for (int i = 0; i < labels.length; i++) {

            int label = labels[i];
            if (count.containsKey(label)) {

                count.put(label, count.get(label) + 1);
            } else {
                count.put(label, 1);

            }

        }

        return count;
    }

    public double[] getProfileX(Atom[] atoms) {
        double[] p = new double[atoms.length];
        for (int i = 0; i < atoms.length; i++) {
            p[i] = atoms[i].getX();
        }
        return p;
    }

    public double[] getProfileZ(Atom[] atoms) {
        double[] p = new double[atoms.length];
        for (int i = 0; i < atoms.length; i++) {
            p[i] = atoms[i].getZ();
        }
        return p;
    }

    public double[] getProfileY(Atom[] atoms) {
        double[] p = new double[atoms.length];
        for (int i = 0; i < atoms.length; i++) {
            p[i] = atoms[i].getY();
        }
        return p;
    }

    public void rotateShiftAtoms(Atom[] ca, Matrix rotMatrix, Atom tranMatrix) {

        for (int i = 0; i < ca.length; i++) {
            // Atom c = (Atom) ca[i].clone();
            Calc.rotate(ca[i], rotMatrix);
            Calc.shift(ca[i], tranMatrix);
            // ca[i] = c;
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

    public Atom genRandomTran() {

        Random random = new Random();
        Atom tranMatrix = atoms[0];
        double rangeMin = 1;
        double rangeMax = 50;

        double ang_x = rangeMin + (rangeMax - rangeMin) * random.nextDouble();
        double ang_y = rangeMin + (rangeMax - rangeMin) * random.nextDouble();
        double ang_z = rangeMin + (rangeMax - rangeMin) * random.nextDouble();

        tranMatrix.setX(ang_x);
        tranMatrix.setY(ang_y);
        tranMatrix.setZ(ang_z);

        return tranMatrix;
    }
}
