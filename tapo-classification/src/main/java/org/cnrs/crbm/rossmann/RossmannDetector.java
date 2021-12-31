package org.cnrs.crbm.rossmann;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.SuperimposerOutput;
import org.cnrs.crbm.lib.repeats.module.Pro2Vector;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.VectorShape;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 8/19/2015.
 */
public class RossmannDetector {


    public static void main(String[] args) throws Exception {


//        RossmannDetector detector = new RossmannDetector("4kbf_A");
//        RossmannDetector detector = new RossmannDetector("1z63_A");
//        RossmannDetector detector = new RossmannDetector("3mwy_W");
//        RossmannDetector detector = new RossmannDetector("1ad3_A");
//        detector.findRossmannRegion();

        RossmannDetector detector = new RossmannDetector("1efk_A");
        detector.findRossmannRegion();


//        try {
//
//            List<String> lines = DataIO.readLines("data/rossmann/rossmannPdbList.in");
//
//            for (String line : lines) {
//                RossmannDetector detector = new RossmannDetector(line);
//                detector.findRossmannRegion();
//            }
//        } catch (StructureException e) {
//            e.printStackTrace();
//        }


    }

    private String pdb;

    public RossmannDetector(String pdb) {
        this.pdb = pdb;
    }

    public void findRossmannRegion() throws Exception {

        String pdbCode = pdb.substring(0, 4);
        String pdbChain = pdb.substring(5, 6);

        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Atom[] atoms = repeatFinder.getAtoms();
        String strSS = repeatFinder.getStrSS();

        Pro2Vector pro2Vector = new Pro2Vector();


//        List<ProVector> vectors = ProVector.toSecondaryVector(repeatFinder.getFeatures()
//                .getVectors());

        List<ProVector> vectors = ProVector.toSecondaryVector(VectorShape
                .getVectors(atoms, strSS));


        int winsize = 4;

        List<List<ProVector>> list = new ArrayList<List<ProVector>>();
        for (int i = 0; i < vectors.size() - winsize + 1; i++) {
            List<ProVector> sub_list = new ArrayList<ProVector>();
            for (int j = 0; j < winsize; j++) {
                ProVector v = vectors.get(i + j);
                sub_list.add(v);
            }
            list.add(sub_list);
        }
        double[] consensus = new double[atoms.length];
        int numberOfRegions = 0;
        Superimposer superimposer = new Superimposer();

        for (int i = 0; i < list.size() - 2; i++) {
            List<ProVector> sub_list_ref = list.get(i);
            List<ProVector> sub_list = list.get(i + 2);
            Atom[] as1 = Fragement.getFragementsofAtoms(atoms, sub_list_ref.get(0)
                    .getPosStart(), sub_list_ref.get(sub_list_ref.size() - 1)
                    .getPosEnd());

            Atom[] as2 = Fragement.getFragementsofAtoms(atoms, sub_list.get(0)
                    .getPosStart(), sub_list.get(sub_list.size() - 1)
                    .getPosEnd());
            int startRegion = sub_list_ref.get(0)
                    .getPosStart();
            int endRegion = sub_list.get(sub_list.size() - 1)
                    .getPosEnd();

            // check pattern
            String pattern_sub_list_ref = pro2Vector.getPatterOfVectors(sub_list_ref);
            String pattern_sub_list = pro2Vector.getPatterOfVectors(sub_list);

            if (pattern_sub_list_ref.equals("HBHB") || pattern_sub_list_ref.equals("BHBH")) {

                SuperimposerOutput superOut = superimposer.compareTwoStructures(as1, as2);
//                double vscore = superimposer.compareListVector(sub_list_ref, sub_list);
                double score = superOut.getTmScore();
                if (pattern_sub_list.equals(pattern_sub_list_ref) && score >= 0.4) {
//                System.out.println(tmScore);
                    numberOfRegions++;
                    for (int j = startRegion; j <= endRegion; j++) {
                        consensus[j] += score;

                    }

                }
            }
            //System.out.println(pro2Vector.getPatterOfVectors(sub_list));

        }

        //normalize
//        if (numberOfRegions > 0) {
//            for (int i = 0; i < consensus.length; i++) {
//                consensus[i] = consensus[i] / numberOfRegions;
//                //consensusLen[i] = consensusLen[i] / lstRep.size();
//            }
//        }
//        for (int i = 0; i < consensus.length; i++) {
//            System.out.print(NumberFormatUtils.format(consensus[i]) + ",");
//        }
//        System.out.println();


        List<RossmanRegion> regions = findRegion(consensus);

        if (regions.size() > 0) {
            StringBuffer buffer = new StringBuffer();
            for (RossmanRegion region : regions) {
                int start = PdbTools.getResSeq(atoms, region.getStart());
                int end = PdbTools.getResSeq(atoms, region.getEnd());
                buffer.append(start + "-" + end + ";");
                System.out.println(region.getStart() + "-" + region.getEnd());

            }
            String units = buffer.toString();
            units = units.substring(0, units.length() - 1);
            System.out.println(pdb + "\tmethod\tclus\t1\t" + "0" + "\t" + "0-0" + "\t" + units + "\t0.0");
        }


    }

    /**
     * display pymolreable output format
     *
     * @param lstRegions
     */
    public void displayPymol(List<RossmanRegion> lstRegions) {
        //


    }

    final double THRESHOLD_CUTOFF = 0.4;

    public List<RossmanRegion> findRegion(double[] concescus) {
        List<RossmanRegion> lstRegions = new ArrayList<RossmanRegion>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < concescus.length; i++) {
            double ch = concescus[i];
            if (ch >= THRESHOLD_CUTOFF) {
                start = i;
                while (i < concescus.length
                        && concescus[i] >= THRESHOLD_CUTOFF) {
                    i++;
                }
                i = end = i - 1;
                // save
                RossmanRegion rc = new RossmanRegion(start, end);
                if (rc.size() > 10)
                    lstRegions.add(rc);
            }


        }
        return lstRegions;
    }

}
