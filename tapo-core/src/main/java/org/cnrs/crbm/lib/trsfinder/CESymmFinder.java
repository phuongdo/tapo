package org.cnrs.crbm.lib.trsfinder;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.cnrs.crbm.lib.multalign.MSAWriter;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 5/26/2015.
 */
public class CESymmFinder extends Finder {
    static Logger logger = LoggerFactory.getLogger(CESymmFinder.class);

    public CESymmFinder(Features features) {
        super(features);
        this.name = "CESymm";
    }

    @Override
    public void findRepeat(Features features) {

        try {


            // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
            // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
            //        System.out.println(structure.toString());


            // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
            // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
            //        System.out.println(structure.toString());

            Structure structure = features.getChain().getStructure().clone();
            Atom[] atoms = StructureTools.getRepresentativeAtomArray(structure.getChainByPDB(features.getPdbChain()));
//            Atom[] ca1 = features.getAtoms();
//            Atom[] atoms = features.getAtoms();
            //Initialize the algorithm
            CeSymm ceSymm = new CeSymm();

            //Choose some parameters
            CESymmParameters params = (CESymmParameters) ceSymm.getParameters();
            params.setRefineMethod(CESymmParameters.RefineMethod.SINGLE);
            params.setOptimization(true);
            params.setMultipleAxes(true);
            params.reset();


            //Run the symmetry analysis - alignment as an output
            MultipleAlignment symmetry = ceSymm.analyze(atoms, params);
//            System.out.println(symmetry.getScore("AvgTM-score"));
            int symmNr = symmetry.size();
            //double symmScore = symmNr + symmetry.getScore("AvgTM-score");
            double symmScore = symmetry.getScore("AvgTM-score");
            this.combineScore.setCeScore(symmScore);
//            System.out.println(symmScore + ":" + symmNr);
            if (symmScore >= 0.4) {
                List<List<Integer>> multiples_tmp = MSAWriter.tapoMsaFormat(symmetry);
                Repeat candiate = new Repeat();
                for (List<Integer> lst : multiples_tmp) {
                    //remove -1
                    List<Integer> lstPos = new ArrayList<Integer>();
                    for (int pos : lst) {
                        if (pos > 0) {
                            lstPos.add(pos);
                        }
                    }
                    if (lstPos.size() > 1)
                        candiate.getRepeats().add(new RepeatContent(lstPos.get(0), lstPos.get(lstPos.size() - 1)));
                }

                List<Repeat> splits = this.splitRepeat(candiate);
                for (Repeat split : splits) {
                    if (split.getRepeats().size() >= 2)
                        this.repeats.add(split);
                    // System.out.println(split);
                }

            }


//            Atom[] ca1 = features.getAtoms();
//            Atom[] ca2 = StructureTools.cloneCAArray(ca1);
//            CeSymm ceSymm = new CeSymm();
//            AFPChain afpChain = ceSymm.pairAlign(ca1, ca2);
//            int symmNr = new SequenceFunctionOrderDetector().calculateOrder(afpChain, ca1);
//            double symmScore = symmNr + afpChain.getTMScore();
//            this.combineScore.setCeScore(symmScore);
//
//            if (symmScore >= 1.4) {// from author
//                // repeats
//                int L = 0;// possible repeats length
//                if (symmScore < 2) {
//                    L = ca1.length / 2;
//                } else {
//
//                    L = ca1.length / symmNr;
//                }
//                if (L > 0) {
//                    // generate TR candidates from CE alignment
////                    MutilAlign mutilAlign = new MutilAlign();
////                    Map<Integer, Integer> pairs = mutilAlign.toAlignedPairs(afpChain, 0, 0);
//
////
////
//
//
////                    double avgL = 0.0;
////                    for (Map.Entry<Integer, Integer> entry : pairs.entrySet()) {
//////                        System.out.println(entry.getKey() + " : "
//////                                + entry.getValue());
////                        alg1[index] = entry.getKey();
////                        alg2[index] = entry.getValue();
//////                        avgL += Math.abs(entry.getValue() - entry.getKey());
////
////                        index++;
////                    }
//
//
//                    int[][][] optAln = afpChain.getOptAln();
//                    int[] blockLen = afpChain.getOptLen();
//                    // find the maximum block which is max length;
//                    int indexBlock = 0;
//                    int maxBlockLen = 0;
//
//
//                    for (int block = 0; block < afpChain.getBlockNum(); block++) {
//                        if (maxBlockLen < blockLen[block]) {
//                            indexBlock = block;
//                            maxBlockLen = blockLen[block];
//                        }
//                    }
//
//                    double avgL = 0.0;
//                    int noAlgRes = 0;
//
//                    int[] alg1max = new int[blockLen[indexBlock]];
//                    int[] alg2max = new int[blockLen[indexBlock]];
//
//                    int index = 0;
//                    for (int i = 0; i < blockLen[indexBlock]; i++) {
//                        int posA1 = optAln[indexBlock][0][i];
//                        int posA2 = optAln[indexBlock][1][i];
//                        avgL += Math.abs(posA1 - posA2);
//                        alg1max[index] = posA1;
//                        alg2max[index] = posA2;
//                        index++;
//                        noAlgRes++;
//
//                    }
//                    if (noAlgRes > 0) {
//                        avgL = avgL / noAlgRes;
//                        if (symmNr == 1)
//                            L = (int) avgL;
//                    }
//                    if (L > 10) {
////                            TMEvaluation tmEvaluation = new TMEvaluation();
////                            FinderOutput finderOut = tmEvaluation.findRepeatsBasedOnTMScore((int) L, 0, atoms.length, features);
////                            List<Repeat> lrepeats = finderOut.getRepeats();
////
////                            if (lrepeats.size() > 0)
////                                this.repeats.addAll(lrepeats);
//
////                        Repeat repeat = new Repeat();
//                        int startRegion = -1;
//                        int endRegion = -1;
//                        double predLen = Math.abs(alg1max[0] - alg2max[0]);
//                        if (0.8 * L < predLen && predLen < L * 1.2) {
//                            if (alg1max[0] < alg2max[0]) {
////                                int previousEnd = -1;
////                                for (int i = 0; i < alg1max.length; i = i + L) {
////                                    int start = alg1max[i];
////                                    if (previousEnd > start)
////                                        start = previousEnd + 1;
////                                    int end = start + L - 1;
////                                    previousEnd = end;
////                                    if (end > atoms.length)
////                                        end = atoms.length - 1;
////                                    if (end - start + 1 > L / 2)
////                                        repeat.getRepeats().add(new RepeatContent(start, end));
////                                }
//
//
//                                startRegion = alg1max[0];
//                                endRegion = alg2max[alg2max.length - 1];
//
//
//                            } else {
//                                startRegion = alg2max[0];
//                                endRegion = alg1max[alg2max.length - 1];
//
//                            }
//
//                            // extract repeat candidates
//                            Repeat candiate = new Repeat();
//                            for (int i = startRegion; i < endRegion; i = i + L) {
//                                int start = i;
//                                int end = i + L - 1;
//                                if (end > endRegion - 1)
//                                    end = endRegion - 1;
//                                RepeatContent ru = new RepeatContent(start, end);
//                                if (ru.size() > (L / 2))
//                                    candiate.getRepeats().add(ru);
//                            }
//                            if (candiate.getRepeats().size() >= 2)
//                                this.repeats.add(candiate);
//
////                        if (repeat.getRepeats().size() >= 2) {
////                            repeat.sortByPosition();
////                            this.repeats.add(repeat);
////                        }
//                        }
//                    }
//
//                }
//            }


        } catch (Exception e) {
//            e.printStackTrace();
            //logger.warn("CE-Symm", e);
            logger.error(e.getMessage() + " with pdb " + features.getPdbCode()
                    + "_" + features.getPdbChain());
        }
    }


    public int findNextIndexOfAlg1(int pos, int[] alg1max, int[] alg2max) {
        int index = -1;
        for (int i = 0; i < alg1max.length; i = i + 1) {
            if (pos <= alg1max[i]) {
                index = i;
            } else {
                break;
            }
        }
        return index;

    }


}
