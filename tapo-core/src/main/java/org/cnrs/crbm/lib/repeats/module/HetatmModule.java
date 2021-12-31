package org.cnrs.crbm.lib.repeats.module;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.multalign.TSim;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.seqalign.RepeatSeq;
import org.cnrs.crbm.lib.seqalign.RepeatUnitSeq;
import org.cnrs.crbm.lib.seqalign.Trust;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.*;

/**
 * Created by pdoviet on 5/22/2015.
 */
public class HetatmModule {


    // CA for calcium ION
    // SF4 and SO4 for sulphate ION

    //defaults
    //String[] atomNames = {"CA", "SF4", "SO4", "CMO", "FE2", "ZN"};
    String[] atomNames = null;
    MutilAlign mutilAlign = new MutilAlign();
    Map<String, Double> vanderwaalsradii = new HashMap<String, Double>();

    public double getMaxHScore() {
        return maxHScore;
    }

    double maxHScore = 0.0;

    public HetatmModule() {
        // load atomNames
        loadLigand();
        loadVanderWaal();
    }

    public static void main(String[] args) throws Exception {
        HetatmModule hetatmModule = new HetatmModule();
        String pdbCode = "4dh2";
        String pdbChain = "B";
//        pdbCode = "4i0x";
//        pdbChain = "B";
//        pdbCode = "1sw8";
//        pdbChain = "A";
//        pdbCode = "1y6w";
//        pdbChain = "A";
//        pdbCode = "2ctt";
//        pdbChain = "A";
        pdbCode = "3dd4";
        pdbChain = "A";
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        Chain chain = structure.getChainByPDB(pdbChain);
        Atom[] atoms = PdbTools.getAtomCAArray(structure.getChainByPDB(pdbChain));
        HetatmModule module = new HetatmModule();
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        Features feature = repeatFinder.getFeatures();
//        String seqAp = feature.getStrSeqAp();
//        System.out.println(seqAp);
//        System.out.println(feature.getStrSeq());
        System.out.println(hetatmModule.findRepeats(feature));
    }


    private void loadLigand() {
        List<String> lines = DataIO.readLines("conf/LigandListFunFOLD.txt");
        List<String> ligands = new ArrayList<String>();
        for (String line : lines) {
            if ((!line.startsWith("#")) && (!line.startsWith("//")) &&
                    (line.length() > 0)) {
                String ligand = line.substring(0, line.indexOf("="));
                ligands.add(ligand);
            }
        }
        atomNames = (String[]) ligands.toArray(new String[ligands.size()]);
    }

    private void loadVanderWaal() {
        List<String> lines = DataIO.readLines("conf/VanderWaalsRadius.cfg");
        for (String line : lines) {
            StringTokenizer strToken = new StringTokenizer(line, ",");
            String atomName = strToken.nextToken().toUpperCase();
            Double size = new Double(strToken.nextToken());
            vanderwaalsradii.put(atomName, size);
        }

    }

    public double getHETAMScore(Features features) {
        Atom[] atoms = features.getAtoms();
        Chain chain = features.getChain();
        double hemScore = 0.0;
        List<List<ContactRegion>> listRegions = this.findingContactRegionWithHETAMGroup(chain, atoms);
        for (List<ContactRegion> regions : listRegions) {
            if (regions.size() >= 2) {
                int[] labels = new int[atoms.length];
                // build pattern
                for (ContactRegion region : regions) {

                    String pattern = region.getContactPattern();
                    for (int i = 0; i < pattern.length(); i++) {
                        labels[i + region.getStart()] = Integer.parseInt(pattern.charAt(i) + "");
                    }

                }
                StringBuffer sBu = new StringBuffer();
                for (int i = 0; i < labels.length; i++) {
                    sBu.append(labels[i]);
                }
                String cpattern = sBu.toString();
                System.out.println(cpattern);
                for (int k = 0; k < regions.size() - 1; k++) {
                    // start region
                    int start = regions.get(k).getStart();
                    // end regions
                    int end = regions.get(k + 1).getEnd();
                    int regionLength = end - start + 1;
                    if (start - regionLength < 0)
                        start = 0;
                    else
                        start = start - regionLength;

                    if (end + regionLength > atoms.length)
                        end = atoms.length - 1;
                    else
                        end = end + regionLength;

                    try {
                        Repeat candidate = this.processBYSuperimpose(atoms, regions.get(k), regions.get(k + 1), cpattern);
                        hemScore = Math.max(candidate.getScore(), hemScore);

                    } catch (Exception e) {
//                        e.printStackTrace();
                    }
                }
            }
        }
        return hemScore;
    }


    /**
     * Finding repeats
     *
     * @param features
     * @return
     */
    public List<Repeat> findRepeats(Features features) {
        Atom[] atoms = features.getAtoms();
        Chain chain = features.getChain();
        List<Repeat> repeats = new ArrayList<Repeat>();
        List<List<ContactRegion>> listRegions = this.findingContactRegionWithHETAMGroup(chain, atoms);
        for (List<ContactRegion> regions : listRegions) {

            ContactRegion.sortByPosition(regions);
            if (regions.size() >= 2) {
                int[] labels = new int[atoms.length];
                // build pattern
                for (ContactRegion region : regions) {

                    String pattern = region.getContactPattern();
                    for (int i = 0; i < pattern.length(); i++) {
                        labels[i + region.getStart()] = Integer.parseInt(pattern.charAt(i) + "");
                    }

                }
                StringBuffer sBu = new StringBuffer();
                for (int i = 0; i < labels.length; i++) {
                    sBu.append(labels[i]);
                }
                String cpattern = sBu.toString();
                //System.out.println(cpattern);
                for (int k = 0; k < regions.size() - 1; k++) {
                    // start region
                    int start = regions.get(k).getStart();
                    // end regions
                    int end = regions.get(k + 1).getEnd();
                    int regionLength = end - start + 1;
                    if (start - regionLength < 0)
                        start = 0;
                    else
                        start = start - regionLength;

                    if (end + regionLength > atoms.length)
                        end = atoms.length - 1;
                    else
                        end = end + regionLength;

                    // extract region

//                    String strRegionAp = seqAp.substring(start, end + 1);
//                    String strRegionCpattern = cpattern.substring(start, end + 1);
//                    //double score = module.processByTRUST(strRegionAp, strRegionCpattern);

                    //double score = module.processBYSuperimpose(atoms, regions.get(k), regions.get(k + 1), cpattern);
                    //System.out.println(start + "\t" + end + "\t" + score);

                    try {
                        Repeat candidate = this.processBYSuperimpose(atoms, regions.get(k), regions.get(k + 1), cpattern);
//                        System.out.println(candidate.getScore());
                        maxHScore = Math.max(candidate.getScore(), maxHScore);
                        if (candidate.getScore() > 0.5) {
                            repeats.add(candidate);
                        }

                    } catch (Exception e) {
                        e.printStackTrace();
                    }


                }


            }

        }
        return repeats;
    }


    public Repeat processBYSuperimpose(Atom[] atoms, ContactRegion seed1, ContactRegion seed2, String cpattern) throws Exception {
        // the best candidate

        Repeat theBest = new Repeat();
        try {
            // check 3 possibilities
            List<Repeat> candidates = new ArrayList<Repeat>();
            int L = seed2.getStart() - seed1.getStart() + 1;
            /**
             * ----XXXX------------XXXX--------------
             *    |_______________||_______________|
             *            L               L
             */
            int end2 = seed2.getStart();
            if (end2 + L > atoms.length)
                end2 = atoms.length - 1;
            else
                end2 = end2 + L;

            int startSet1 = seed1.getStart();
            int startSet2 = seed2.getStart();
            Repeat candiate1 = new Repeat();
            candiate1.getRepeats().add(new RepeatContent(startSet1, startSet2 - 1));
            candiate1.getRepeats().add(new RepeatContent(startSet2, end2));

//        Atom[] setAtoms1 = Fragement.getFragementsofAtoms(atoms, startSet1, startSet2 - 1);
//        Atom[] setAtoms2 = Fragement.getFragementsofAtoms(atoms, startSet2, end2);
//        AFPChain afpChain = mutilAlign.pairAlign(setAtoms1, setAtoms2);
//        double tmScore = afpChain.getTMScore();
//        double hScore = this.getScoreH(afpChain, cpattern, startSet1, startSet2);
//        double finalScore = (tmScore + hScore) / 2;
            //System.out.println(">>"+startSet1+"to"+end2);


            /**
             * Middle of the seed
             * ----XXXX------------XXXX--------------
             *|_____________||_______________|
             *      L               L
             */

            /**
             * bug fixed from version 1.1.3
             */
            int start1 = seed1.getStart();
            if (start1 - L / 2 < 0)
                start1 = 0;
            else
                start1 = start1 - L / 2;
            end2 = seed2.getStart();
            if (end2 + L / 2 > atoms.length)
                end2 = atoms.length - 1;
            else
                end2 = end2 + L / 2;

            int mid = (seed1.getEnd() + seed2.getStart()) / 2;

            startSet1 = start1;
            startSet2 = mid;

            Repeat candiate2 = new Repeat();
            candiate2.getRepeats().add(new RepeatContent(startSet1, mid - 1));
            candiate2.getRepeats().add(new RepeatContent(startSet2, end2));

//        setAtoms1 = Fragement.getFragementsofAtoms(atoms, startSet1,mid-1);
//        setAtoms2 = Fragement.getFragementsofAtoms(atoms, startSet2, end2);
//        afpChain = mutilAlign.pairAlign(setAtoms1, setAtoms2);
//        tmScore = afpChain.getTMScore();
//        hScore = this.getScoreH(afpChain, cpattern, startSet1, startSet2);
//        finalScore = Math.max(finalScore, (tmScore + hScore) / 2);


            /**
             *----------XXXX------------XXXX-------------
             *|_____________||______________|
             *     L                L
             */
            start1 = seed1.getStart();
            if (start1 - L < 0)
                start1 = 0;
            else
                start1 = start1 - L;

            startSet1 = start1;
            startSet2 = seed1.getEnd() + 1;

            Repeat candiate3 = new Repeat();
            candiate3.getRepeats().add(new RepeatContent(startSet1, seed1.getEnd()));
            candiate3.getRepeats().add(new RepeatContent(startSet2, seed2.getEnd()));

            candidates.add(candiate1);
            candidates.add(candiate2);
            candidates.add(candiate3);


            double maxScore = 0.0;
            for (Repeat candiate : candidates) {
                RepeatContent ru1 = candiate.getRepeats().get(0);
                RepeatContent ru2 = candiate.getRepeats().get(1);
                Atom[] setAtoms1 = Fragement.getFragementsofAtoms(atoms, ru1.getStart(), ru1.getEnd());
                Atom[] setAtoms2 = Fragement.getFragementsofAtoms(atoms, ru2.getStart(), ru2.getEnd());
                List<Atom> listIncludedHEATM1 = new ArrayList<Atom>();
                List<Atom> listIncludedHEATM2 = new ArrayList<Atom>();

                for (Atom atom : setAtoms1) {
                    listIncludedHEATM1.add(atom);
                }
                for (Atom atom : seed1.getHematom().getAtoms()) {
//                Atom a = new AtomImpl();
//                a.setCoords(atom.getCoords());
//                a.setName(" CA ");
//                a.setGroup((Group) atom.getGroup().clone());

                    Atom a = (Atom) atoms[0].clone();
                    a.setCoords(atom.getCoords());
                    a.setGroup((Group) atoms[0].getGroup().clone());
                    listIncludedHEATM1.add(a);
                }

                for (Atom atom : setAtoms2) {
                    listIncludedHEATM2.add(atom);
                }
                for (Atom atom : seed2.getHematom().getAtoms()) {
                    Atom a = (Atom) atoms[0].clone();
                    a.setGroup((Group) atoms[0].getGroup().clone());
                    a.setCoords(atom.getCoords());
                    listIncludedHEATM2.add(a);
                }

                setAtoms1 = (Atom[]) listIncludedHEATM1.toArray(new Atom[listIncludedHEATM1.size()]);
                setAtoms2 = (Atom[]) listIncludedHEATM2.toArray(new Atom[listIncludedHEATM2.size()]);
                AFPChain afpChain = mutilAlign.pairAlign(setAtoms1, setAtoms2);
                double tmScore = afpChain.getTMScore();
                double hScore = this.getScoreH(afpChain, cpattern, ru1.getStart(), ru2.getStart());
                double finalScore = (tmScore + hScore) / 2;
                if (finalScore > maxScore) {
                    maxScore = finalScore;
                    theBest = candiate;
                    theBest.setScore(maxScore);
                }

            }


//        setAtoms1 = Fragement.getFragementsofAtoms(atoms, startSet1,seed1.getEnd());
//        setAtoms2 = Fragement.getFragementsofAtoms(atoms, startSet2, seed2.getEnd());
//        afpChain = mutilAlign.pairAlign(setAtoms1, setAtoms2);
//        tmScore = afpChain.getTMScore();
//        hScore = this.getScoreH(afpChain, cpattern, startSet1, startSet2);
//        finalScore = Math.max(finalScore, (tmScore + hScore) / 2);
        } catch (Exception ex) {
            //ex.printStackTrace();
        }
        return theBest;
    }


    private double getScoreH(AFPChain afpChain, String cpattern, int startSet1, int startSet2) {

        Map<Integer, Integer> alignPairs = mutilAlign.toAlignedPairs(afpChain,
                startSet1, startSet2);

        int leng = alignPairs.size();
        int totalH = 0;
        int totalPairsH = 0;
        for (Map.Entry<Integer, Integer> entry : alignPairs.entrySet()) {
            char ch1 = cpattern.charAt(entry.getKey());
            char ch2 = cpattern.charAt(entry.getValue());
            if (ch1 == '1' || ch2 == '1') {
                totalH++;

                if (ch1 == ch2) {
                    totalPairsH++;
                }

            }
        }

        if (totalH > 0) {
            return (double) totalPairsH / totalH;
        }
        return 0.0;
    }


    public double processByTRUST(String strRegionAp, String strRegionCpattern) {

        // process by TRUST
        int gapo = -6;
        int gape = -2;
        Trust trust = new Trust("conf/BLOSUM_COMBINED", "EKLIBVPSDCGHAFXYOWMN",
                "", strRegionAp, gapo, gape);

//        Trust trust = new Trust("conf/BLOSUM62", "ABCDEFGHIKLMNPQRSTVWXYZO",
//                "", strRegionAp, gapo, gape);
        //        Collection repeats_trust = trust.getRepeatsOutput();
        //List<String> msa = trust.getMsa();
        double finalScore = 0.0;
        try {

            List<RepeatSeq> repeats = trust.getRepeats();
            for (RepeatSeq repeat : repeats) {
                List<List<Integer>> multiples = new ArrayList<List<Integer>>();
                int msaLeng = repeat.getUnits().get(0).size();
                for (RepeatUnitSeq ru : repeat.getUnits()) {
                    List<Integer> list = new ArrayList<Integer>();
                    int start = ru.getStart();
                    String seqs = ru.getStrAlign();
                    int k = 0;
                    for (int j = 0; j < msaLeng; j++) {
                        char AA = seqs.charAt(j);
                        int pos = -1;
                        if (AA != '-') {
                            pos = start + k;
                            k++;
                        } else {
                            pos = -1;
                        }
                        list.add(pos);
                    }
                    multiples.add(list);
                    System.out.println(ru.getStrAlign() + "\t" + ru.getStart() + "-" + ru.getEnd());
                }// END FOR

                // process multiple
                List<String> cpatternMSA = mutilAlign.convertToMSA(multiples, strRegionCpattern);
                List<String> msa = mutilAlign.convertToMSA(multiples, strRegionAp);

                TSim tsim = new TSim();
//                System.out.println(tsim.getTsimScore(msa));
//                System.out.println(tsim.getHETATMScore(cpatternMSA));
                finalScore = (tsim.getTsimScore(msa) + tsim.getHETATMScore(cpatternMSA)) / 2;


            }
        } catch (Exception ex) {

            ex.printStackTrace();
            // do nothing here
        }

        return finalScore;

    }


    public List<List<ContactRegion>> findingContactRegionWithHETAMGroup(Chain chain, Atom[] caAtoms) {
        // get hetatm from one chain
        List<List<Group>> listHETATMs = new ArrayList<List<Group>>();

        List<Group> groups = chain.getAtomGroups(GroupType.HETATM);
        for (int i = 0; i < atomNames.length; i++) {
            // a temp container for the atoms of this group
            //List<Atom> thisGroupAtoms = new ArrayList<Atom>();
            List<Group> molecularheatm = new ArrayList<Group>();
            String atomName = atomNames[i];

            for (Group g : groups) {
                // flag to check if this group contains all the requested atoms.
                //boolean thisGroupAllAtoms = true;
                //System.out.println(g.getPDBName());
                if (g.getPDBName().trim().equals(atomName)) {
                    molecularheatm.add(g);
                }

            }

            if (molecularheatm.size() > 0)
                listHETATMs.add(molecularheatm);
        }
        List<List<ContactRegion>> listRegions = new ArrayList<List<ContactRegion>>();
        List<ContactRegion> regions = new ArrayList<ContactRegion>();
        for (List<Group> thisGroupAtoms : listHETATMs) {

            // for one type of HEATM

            Map<Integer, Group> mapheatmPosition = new HashMap<Integer, Group>();

            for (Group molecularH : thisGroupAtoms) {
                int[] labels = new int[caAtoms.length];
                int startRegion = 0;
                int endRegion = 0;
                for (int i = 0; i < caAtoms.length; i++) {
                    Group residue = caAtoms[i].getGroup();
                    try {

                        for (Atom atom : residue.getAtoms()) {
                            for (Atom atomH : molecularH.getAtoms()) {

                                double dist = Calc.getDistance(atom, atomH);
                                String nameAtom = atom.getName().trim();
                                String nameAtomH = atomH.getName().trim();
                                double sizeAtom = 0.0;
                                double sizeAtomH = 0.0;
                                if (vanderwaalsradii.containsKey(nameAtom)) {
                                    sizeAtom = vanderwaalsradii.get(nameAtom);
                                }
                                if (vanderwaalsradii.containsKey(nameAtomH)) {
                                    sizeAtomH = vanderwaalsradii.get(nameAtomH);
                                }

                                if (dist < sizeAtom + sizeAtomH + 0.5D) {

                                    labels[i] = 1;
                                    endRegion = Math.max(i, endRegion);
                                    mapheatmPosition.put(i, molecularH);
                                    break;
                                }

                            }
                        }

                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                }

                for (int i = 0; i < labels.length; i++) {
                    if (labels[i] == 1) {
                        startRegion = i;
                        break;
                    }

                }
                // extract region nows

                /**
                 * we consider the contact region there are at least 2 residues in contacting.
                 */
                if (startRegion > 0 & endRegion >= startRegion) {
                    ContactRegion contactRegion = new ContactRegion();
                    StringBuffer buffer = new StringBuffer();
                    int count = 0;
                    for (int k = startRegion; k <= endRegion; k++) {
                        buffer.append(labels[k]);
                        count++;
                    }
                    if (count >= 2) {
                        contactRegion.setHematom(molecularH);
                        contactRegion.setStart(startRegion);
                        contactRegion.setEnd(endRegion);
                        contactRegion.setContactPattern(buffer.toString());
                        regions.add(contactRegion);
                    }

                }


            }
            if (regions.size() > 0)
                listRegions.add(regions);


//            System.out.println();

        }// end group
        return listRegions;
    }// END FOR


    public List<ContactRegion> etractRegion(int[] labels, Map<Integer, Group> mapheatmPosition) {

        List<ContactRegion> regions = new ArrayList<ContactRegion>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < labels.length; i++) {
            int ch = labels[i];
            if (ch == 1) {
                start = i;
                while (i < labels.length
                        && labels[i] == 1) {
                    i++;
                }
                i = end = i - 1;
                // save

                if (end - start + 1 >= 3) {
                    //System.out.println(start + "-" + end);
                    ContactRegion contactRegion = new ContactRegion();
                    contactRegion.setStart(start);
                    contactRegion.setEnd(end);
                    regions.add(contactRegion);
                    // set HEATM to a region
                    for (int j = start; j <= end; j++) {
                        if (mapheatmPosition.containsKey(j)) {
                            contactRegion.setHematom(mapheatmPosition.get(j));
                            break;
                        }

                    }
                }

            }

        }// END FOR

        return regions;

    }


}

class ContactRegion {


    public static void sortByPosition(List<ContactRegion> regions) {
        Collections.sort(regions, new Comparator<ContactRegion>() {
            public int compare(ContactRegion s1, ContactRegion s2) {
                if (s1.getStart() > s2.getStart())
                    return 1;
                else
                    return -1;
            }
        });

    }

    private int start;

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    String contactPattern = "";

    public String getContactPattern() {
        return contactPattern;
    }

    public void setContactPattern(String contactPattern) {
        this.contactPattern = contactPattern;
    }

    private int end;
    private Group hematom;

    public Group getHematom() {
        return hematom;
    }

    public void setHematom(Group hematom) {
        this.hematom = hematom;
    }

    @Override
    public String toString() {
        return this.start + "-" + this.end;
    }
}
