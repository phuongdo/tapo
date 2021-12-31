//package org.cnrs.crbm.lib.trsfinder;
//
//        import org.apache.log4j.Logger;
//
//        import org.biojava.nbio.structure.Atom;
//        import org.biojava.nbio.structure.StructureException;
//        import org.cnrs.crbm.lib.math.PeakDetector;
//        import org.cnrs.crbm.lib.repeats.Fragement;
//        import org.cnrs.crbm.lib.repeats.Superimposer;
//        import org.cnrs.crbm.lib.repeats.module.VectorShape;
//        import org.cnrs.crbm.lib.utils.NumberFormatUtils;
//
//        import java.util.ArrayList;
//        import java.util.List;
//
///**
// * THis class was designed to deal with two peak
// *
// * @author pdoviet
// */
//@Deprecated
//public class RMSDFinder extends Finder {
//
//    final static boolean debug = false;
//
//    final static int MIN_WINSIZE = 20;
//    final static int MAX_WINSIZE = 50;
//    static Logger logger = Logger.getLogger(RMSDFinder.class);
//
//    public RMSDFinder(Features features) {
//        super(features);
//        this.name = "RMSD";
//
//    }
//
//    @Override
//    public void findRepeat(Features features) {
//        Superimposer superimposer = new Superimposer();
//        Atom[] atoms = features.getAtoms();
//        String strSS = features.getStrSS();
//        // winsize = 30, 40, 50;
//        // int winsize = 30;
//        for (int winsize = MIN_WINSIZE; winsize < MAX_WINSIZE; winsize = winsize + 10) {
//            try {
//                List<List<Integer>> seeds = this.scanSeed(features, winsize);
//                for (List<Integer> seed : seeds) {
//                    if (seed.size() == 2) {
//                        int p0 = seed.get(0);
//                        int p1 = seed.get(1);
//                        int p2 = p1 + p1 - p0;
//                        if (p2 >= atoms.length)
//                            p2 = atoms.length - 1;
//                        Atom[] as1 = Fragement.getFragementsofAtoms(atoms, p0,
//                                p1);
//                        Atom[] as2 = Fragement.getFragementsofAtoms(atoms, p1,
//                                p2);
//                        String pattern1 = VectorShape.getSSPattern(strSS
//                                .substring(p0, p1));
//                        String pattern2 = VectorShape.getSSPattern(strSS
//                                .substring(p1, p2));
//
//                        // System.out.println(p0 + "-" + p1 + "-" + p2);
//
//                        // check both structure and pattern
//
//                        boolean check = true;
//                        if (as1.length < 40)
//                            check = pattern1.equals(pattern2);
//
//                        if (superimposer.isTheSameStructure(as1, as2)
//                                && pattern1.length() >= 2
//                                && pattern2.length() >= 2 && check) {
//                            // System.out.println(pattern1 + ":" + pattern2);
//                            org.cnrs.crbm.lib.trsfinder.Repeat aRepeat = new org.cnrs.crbm.lib.trsfinder.Repeat();
//                            aRepeat.getRepeats().add(
//                                    new RepeatContent(p0, p1 - 1));
//                            aRepeat.getRepeats().add(new RepeatContent(p1, p2));
//                            // System.out.println(p0 + "-" + p1 + "-" + p2);
//                            repeats.add(aRepeat);
//
//                        }
//
//                        // other cases
//                        int leng = p1 - p0;
//                        int start1 = p0 - (int) leng / 2;
//                        if (start1 < 0)
//                            start1 = 0;
//                        int end1 = (int) (p1 + p0) / 2;
//                        int start2 = end1 + 1;
//                        int end2 = start2 + leng - 1;
//                        if (end2 >= atoms.length)
//                            end2 = atoms.length - 1;
//
//                        as1 = Fragement.getFragementsofAtoms(atoms, start1,
//                                end1);
//                        as2 = Fragement.getFragementsofAtoms(atoms, start2,
//                                end2);
//                        pattern1 = VectorShape.getSSPattern(strSS.substring(
//                                start1, end1));
//                        pattern2 = VectorShape.getSSPattern(strSS.substring(
//                                start2, end2));
//
//                        // System.out.println(p0 + "-" + p1 + "-" + p2);
//
//                        check = true;
//                        if (as1.length < 40)
//                            check = pattern1.equals(pattern2);
//
//                        // check both structure and pattern
//                        if (superimposer.isTheSameStructure(as1, as2)
//                                && pattern1.length() >= 2
//                                && pattern2.length() >= 2 && check) {
//                            // System.out.println(pattern1 + ":" + pattern2);
//                            org.cnrs.crbm.lib.trsfinder.Repeat aRepeat = new org.cnrs.crbm.lib.trsfinder.Repeat();
//                            aRepeat.getRepeats().add(
//                                    new RepeatContent(start1, end1));
//                            aRepeat.getRepeats().add(
//                                    new RepeatContent(start2, end2));
//                            // System.out.println(p0 + "-" + p1 + "-" + p2);
//                            repeats.add(aRepeat);
//
//                        }
//
//                    } else if (seed.size() > 2) {
//
//                        int p0 = seed.get(0);
//                        int p1 = seed.get(1);
//                        org.cnrs.crbm.lib.trsfinder.Repeat aRepeat = new org.cnrs.crbm.lib.trsfinder.Repeat();
//                        // aRepeat.getRepeats().add(new RepeatContent(p0, p1 -
//                        // 1));
//                        Atom[] as1 = Fragement.getFragementsofAtoms(atoms, p0,
//                                p1);
//                        String pattern1 = VectorShape.getSSPattern(strSS
//                                .substring(p0, p1));
//
//                        for (int i = 1; i < seed.size() - 1; i++) {
//                            int p_i = seed.get(i) + 1;
//                            int p_i_next = seed.get(i + 1);
//
//                            Atom[] as2 = Fragement.getFragementsofAtoms(atoms,
//                                    p_i, p_i_next);
//                            String pattern2 = VectorShape.getSSPattern(strSS
//                                    .substring(p_i, p_i_next));
//                            // check both structure and pattern
//                            if (pattern2.length() >= 2
//                                    && pattern1.length() >= 2
//                                    && superimposer
//                                    .isTheSameStructure(as1, as2)
//
//                                    ) {
//
//                                if (i == 1) {
//                                    aRepeat.getRepeats().add(
//                                            new RepeatContent(p0, p1 - 1));
//                                }
//                                aRepeat.getRepeats().add(
//                                        new RepeatContent(p_i, p_i_next));
//                                // System.out.println(p0 + "-" + p1 + "-" + p2);
//                            } else {
//
//                                // we only consider to the TRs
//                                break;
//                            }
//                        }
//                        if (aRepeat.getRepeats().size() >= 2)
//                            repeats.add(aRepeat);
//
//                    }
//
//                }
//
//            } catch (Exception e) {
//                // TODO Auto-generated catch block
//                logger.error(e.getMessage() + " with pdb "
//                        + features.getPdbCode() + "_" + features.getPdbChain());
//
//            }
//        }
//
//        // add 11/08/2014 for testing
//        //repeats.addAll(this.getPureRMSDFinder(features));
//    }
//
//    public List<Repeat> getPureRMSDFinder(Features features) {
//
//        List<Repeat> pures = new ArrayList<Repeat>();
//        Atom[] atoms = features.getAtoms();
//        // segmentation
//
//        if (atoms.length < 500) {
//            int[] wlist = new int[]{200};
//            for (Integer ssize : wlist) {
//                for (int i = 0; i < atoms.length; ) {
//                    try {
//                        int start = i;
//                        int end = i + ssize;
//
//                        if (end > atoms.length)
//                            end = atoms.length;
//
//                        if ((double) (end - start + 1) / ssize < 0.2)
//                            break;
//                        //System.out.println("X");
//                        for (int winsize = 25; winsize < 100; winsize = winsize + 3) {
//                            //System.out.println(winsize);
////                            if (winsize == 50)
////                                System.out.println();
//                            List<Repeat> lrepeat = AtomFinder.getPureRMSDFinder(features, winsize, start, end);
//                            if (lrepeat.size() > 0) pures.addAll(lrepeat);
//                        }
//
//                        //System.out.println("Y");
//                        i = i + ssize;
//                    } catch (Exception ex) {
//                        ///
//                    }
//
//                }
//            }
//        }
//
//        return pures;
//
//    }
//
//    public List<List<Integer>> scanSeed(Features features, int winsize)
//            throws Exception {
//
//        if (debug)
//            System.out.println("****" + winsize + "********");
//        Superimposer superimposer = new Superimposer();
//        List<List<Integer>> seeds = new ArrayList<List<Integer>>();
//        Atom[] atoms = features.getAtoms();
//        // get seed fragment
//
//        for (int i = 0; i < atoms.length - winsize; i = i + winsize) {
//            Atom[] seedAtoms = Fragement.getFragementsofAtoms(atoms, i, i
//                    + winsize);
//            List<Double> signals = new ArrayList<Double>();
//            for (int j = i + winsize; j < atoms.length - winsize; j++) {
//                Atom[] compareAtoms = Fragement.getFragementsofAtoms(atoms, j,
//                        j + winsize);
//
//                double rmsd = 10;
//                try {
//                    rmsd = superimposer.superimposeSimple(seedAtoms,
//                            compareAtoms);
//                    signals.add(rmsd);
//                    if (debug)
//                        System.out.print(NumberFormatUtils.format(rmsd) + ", ");
//                } catch (StructureException e) {
//                    e.printStackTrace();
//                }
//
//            }
//            if (debug)
//                System.out.println();
//
//            List<Integer> seedPos = new ArrayList<Integer>();
//            seedPos.add(i);
//            if (signals.size() > 0) {
//                int[] downLocations = PeakDetector.getBottoms(signals);
//
//                for (int k = 0; k < downLocations.length; k++) {
//                    // System.out.println(downLocations[k]);
//                    if (signals.get(downLocations[k]) < 6.0)
//                        seedPos.add(i + winsize + downLocations[k]);
//                }
//            }
//
//            // DEBUG>.....
//            if (debug) {
//                System.out.println("****seed********");
//                for (int pos : seedPos) {
//                    System.out.print(pos + " ");
//
//                }
//                System.out.println("");
//            }
//            if (seedPos.size() >= 2) {
//                seeds.add(seedPos);
//
//            }
//
//        }
//
//        if (debug)
//            System.out.println();
//        return seeds;
//
//    }
//}
