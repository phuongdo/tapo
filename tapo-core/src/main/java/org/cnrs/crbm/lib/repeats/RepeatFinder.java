package org.cnrs.crbm.lib.repeats;

import com.google.common.collect.Lists;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.cdhit.CDHIT;
import org.cnrs.crbm.lib.cm.ContactMap;
import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.dssp.DSSP;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.repeats.clusters.ClusterRepeat;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.VectorModule;
import org.cnrs.crbm.lib.sadb.*;
import org.cnrs.crbm.lib.trsfinder.*;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.ml.SVMWapper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class RepeatFinder {
    /**
     * properties
     */
    String pdbCode;
    String pdbChain;
    String output;
    String strSS;
    String strAcc;
    String strSeqAp;
    String strSeqSADB;
    String strSeq;
    //String cmPattern;
    Atom[] atoms;
    List<ProVector> vectors = new ArrayList<ProVector>();
    double[] cmHistos;
    String mode = "";


    Features features = null;

    boolean isTRs = false;

    StringBuilder contentBuilder = new StringBuilder();
    StringBuilder combinedBuilder = new StringBuilder();
    CombineScore combineScore = new CombineScore();

    double svmScore = 0.0;
    boolean tprefined = true;


    static Logger logger = LoggerFactory.getLogger(RepeatFinder.class);

    public Map<Integer, List<Integer>> getCmapsList() {
        return cmapsList;
    }

    private Map<Integer, List<Integer>> cmapsList = new HashMap<Integer, List<Integer>>();

    /**
     * threading support
     */

    private int nrOfProcessors = 1;


    //private final static double QA_THRES = 0.17;

    /**
     * debug
     */
    static boolean DEBUG_PRINT_TRACES = false;


    /**
     * output absolute score
     */

    private boolean outputRelativeScore = true;


    /**
     * Main for TESTING
     *
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {

        // only protein code was provided
        // 1ojvB
        // 2KIV
        // 3lyc
        // 2xgf
        // 3iox

        String pdbCode = "2y7c";
//        String pdbCode = "1oe4";
        String pdbChain = "A";
        if (args.length > 0) {
            pdbCode = args[0];
            pdbChain = args[1];
        }
        findProteinWithChain(pdbCode, pdbChain);
        //findProteinWithCode(pdbCode);

    }

    public static void findProteinWithChain(String pdbCode, String pdbChain) {
        System.out.println("Analyse protein : " + pdbCode + " with chain: "
                + pdbChain);

        try {
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            repeatFinder.show();
            //repeatFinder.findRepeat();
//            System.out.println(repeatFinder.getOutput());

        } catch (Exception ex) {

            //ex.printStackTrace();
            logger.error(ex.getMessage() + " with pdb " + pdbCode + "_"
                    + pdbChain);
        }

    }

    public static void findProteinWithCode(String pdbCode) {

        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        CDHIT cdhit = new CDHIT(structure, 0.9, 5);
        List<String> chains = cdhit.getUniqueChain();
        System.out.println("Analyse protein : " + pdbCode);
        for (String pdbChain : chains) {
            System.out.println(pdbChain);
            try {
                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);

                // repeatFinder.show();
                repeatFinder.findRepeat();
                // System.out.println(repeatFinder.getOutput());
            } catch (Exception ex) {

                ex.printStackTrace();
            }
        }
    }

    public RepeatFinder(String pdbCode, String pdbChain) {

        this.pdbCode = pdbCode;
        this.pdbChain = pdbChain;

        try {
            this.initial();
        } catch (Exception ex) {
            // ex.printStackTrace();
        }
    }

    public RepeatFinder(String pdbCode, String pdbChain, int nCore) {

        this.pdbCode = pdbCode;
        this.pdbChain = pdbChain;
        this.nrOfProcessors = nCore;

        try {
            this.initial();
        } catch (Exception ex) {
             ex.printStackTrace();
        }
    }

    /**
     * Finding tandem repeats in file with PDB format
     *
     * @param pdbfileDir
     * @param nCore
     */
    public RepeatFinder(String pdbfileDir, String pdbCode, String pdbChain, int nCore) {

        this.pdbCode = pdbCode;
        this.pdbChain = pdbChain;
        this.nrOfProcessors = nCore;

        try {
            this.initial(pdbfileDir);
        } catch (Exception ex) {
            // ex.printStackTrace();
        }
    }

    private void initial(String pdbfileDir) throws StructureException {

        // get information Atoms
        // get information of DSSP
        // version 1.0.2 . support get structure directly from local install of
        // PDB database.

        Structure structure = PdbTools.getStructureFromFile(pdbfileDir);

        // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
        // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
        //
        Chain chain = structure.getChainByPDB(pdbChain);
        atoms = PdbTools.getAtomCAArray(chain);

        if (atoms.length < 10)
            throw new StructureException("The number of atoms is too small!!!");
        else if (atoms.length > 10000)

            throw new StructureException("The number of atoms is to big!!!");

        DSSP dssp = new DSSP(pdbfileDir, pdbCode, pdbChain);

        Sequence3D sequence3D = new Sequence3D();

        // secondary structure
        // strSS = dssp.filterSS(dssp.getSS(atoms));
        strSS = dssp.filterSS(dssp.getCombinedSS(atoms, this.pdbChain));

        // sequence
        // strSeq = dssp.getAminoAcidSequence(pdbChain);
//        strSeq = chain.getSeqResSequence();
//        System.out.println(PdbTools.getSeqResSequence(atoms));
        strSeq = PdbTools.getSeqResSequence(atoms);
        // accessible surface
        strAcc = dssp.getRelativeACC(pdbChain, atoms);

        StrucState[] strucStates = dssp.getStrucState(pdbChain, atoms);
        // alphabet
        strSeqAp = sequence3D.getSA(strucStates, strAcc);

        Atom[] cbAtoms = PdbTools.getAtomCBArray(atoms);
        ContactMap contactMap = new ContactMap(cbAtoms, 7.0, dssp.getSS(pdbChain,
                atoms));

        // contact maps
        cmHistos = contactMap.getHistogram();
        cmapsList = contactMap.getCmapsList();

        StructAlphabet sadb = new StructAlphabet();
        // get vectors

        strSeqSADB = sadb.getSA(strucStates);

        VectorModule vectorModule = new VectorModule();
        //vectors = VectorShape.getVectors(atoms, strSS);
        /**
         * MAJOR UPDATE from 26 May 2015
         */
        vectors = vectorModule.extractVectors(atoms, strSS);
        /**
         * END
         */
        features = new Features(pdbCode, pdbChain, "", strSS, strAcc, strSeqAp,
                strSeq, atoms, vectors, cmHistos);
        features.setStrSeqSADB(strSeqSADB);
        features.setChain(chain);
        features.setCmapsList(cmapsList);
//        repeatExtend = new RepeatExtend(features);
    }


    private void initial() throws StructureException {

        // get information Atoms
        // get information of DSSP
        // version 1.0.2 . support get structure directly from local install of
        // PDB database.
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);

        // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
        // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
        //
        Chain chain = structure.getChainByPDB(pdbChain);
        atoms = PdbTools.getAtomCAArray(chain);

        if (atoms.length < 10)
            throw new StructureException("The number of atoms is too small!!!");
        else if (atoms.length > 5000)
            throw new StructureException("The number of atoms is to big!!!");

        DSSP dssp = new DSSP(pdbCode, pdbChain);

        Sequence3D sequence3D = new Sequence3D();

        // secondary structure
        // strSS = dssp.filterSS(dssp.getSS(atoms));
        strSS = dssp.filterSS(dssp.getCombinedSS(atoms, this.pdbChain));

        // sequence
        // strSeq = dssp.getAminoAcidSequence(pdbChain);
//        strSeq = chain.getSeqResSequence();
//        System.out.println(PdbTools.getSeqResSequence(atoms));
        strSeq = PdbTools.getSeqResSequence(atoms);
        // accessible surface
        strAcc = dssp.getRelativeACC(pdbChain, atoms);

        StrucState[] strucStates = dssp.getStrucState(pdbChain, atoms);
        // alphabet
        strSeqAp = sequence3D.getSA(strucStates, strAcc);

        Atom[] cbAtoms = PdbTools.getAtomCBArray(atoms);
        ContactMap contactMap = new ContactMap(cbAtoms, 7.0, dssp.getSS(pdbChain,
                atoms));

        // contact maps
        cmHistos = contactMap.getHistogram();
        cmapsList = contactMap.getCmapsList();
//        for(int i = 0 ; i < cmHistos.length; i ++){
//            System.out.print(cmHistos[i]+",");
//        }
//        System.out.println();

//        cmHistos = contactMap.getHistogramWithRange((int) atoms.length );

//        StringBuffer buffer = new StringBuffer();
//        //build cm pattern
//        for (int i = 0; i < cmHistos.length; i++) {
//            if (cmHistos[i] < (0.3 * cmHistos.length)) {
//                buffer.append(0);
//            } else {
//                buffer.append(1);
//            }
//
//        }
//
//        cmPattern = buffer.toString();

        StructAlphabet sadb = new StructAlphabet();
        // get vectors

        strSeqSADB = sadb.getSA(strucStates);

        VectorModule vectorModule = new VectorModule();
        //vectors = VectorShape.getVectors(atoms, strSS);
        /**
         * MAJOR UPDATE from 26 May 2015
         */
        vectors = vectorModule.extractVectors(atoms, strSS);
        /**
         * END
         */
        features = new Features(pdbCode, pdbChain, "", strSS, strAcc, strSeqAp,
                strSeq, atoms, vectors, cmHistos);
        features.setStrSeqSADB(strSeqSADB);
        features.setChain(chain);
        features.setCmapsList(cmapsList);
//        repeatExtend = new RepeatExtend(features);
    }

    public void findByOneMethod() {
        Finder finder = new TReksFinder(features);
        List<Repeat> repeats = finder.getRepeats();

        if (finder.getNrOfTRs() >= 2)
            System.out.println(pdbCode + "\tyes");
        else
            System.out.println(pdbCode + "\tno");
    }

    public void show() throws Exception {
        Sequence3D sequence3D = new Sequence3D();

        System.out.println(strSeq);
        System.out.println(strSeqAp);
        Seq3d seq3D = new Seq3d(sequence3D.getSAV2NoInOutSide(strSeqAp));
        System.out.println(sequence3D.getSAV2NoInOutSide(strSeqAp));
        System.out.println(seq3D.getSeq());
        System.out.println(strAcc);
        System.out.println(strSeqSADB);
        System.out.println(strSS);

        RMSDSubFinder finder = new RMSDSubFinder(features);

        List<Integer> listPosL = new ArrayList<Integer>();
//        for (int i = 25; i < 60; i = i + 1) {
//            listPosL.add(i);
//        }

//        listPosL.add(33);
        finder.setListWinsize(listPosL);
        System.out.println("start finding");
        finder.start();
        System.out.println("end finding");

//        System.out.println("start extending");
//
//        List<Repeat> repeats = finder.getRepeats();
//        for (Repeat repeat : repeats) {
//            // System.out.println(repeat);
//            try {
//
//                List<RepeatContent> templateTRs = new ArrayList<RepeatContent>();
//                for (RepeatContent r : repeat.getRepeats()) {
//                    templateTRs.add(r);
//                }
//                repeatExtend.extendRepeatRight(repeat, templateTRs);
//                repeatExtend.extendRepeatLeft(repeat, templateTRs);
//
//                System.out.println("***");
//                System.out.println(repeat);
//
//            } catch (Exception ex) {
//                ex.printStackTrace();
//            }
//        }
//        System.out.println("end extending");


        System.out.println(strSS);
        List<Repeat> repeats = finder.getRepeats();
        for (Repeat repeat : repeats) {
            System.out.println(this.showRepeatsByConsole(repeat));
        }
        System.out
                .println(finder.isTRs() + " and nrTRs:" + finder.getNrOfTRs());
        // if(ra)

    }

    public String showRepeatsByConsole(Repeat repeat) {

        StringBuffer buffer = new StringBuffer();
        int[] lables = new int[atoms.length];


        int index = 1;
        for (RepeatContent rc : repeat.getRepeats()) {
            for (int i = rc.getStart(); i <= rc.getEnd(); i++) {
                lables[i] = index;
            }
            index++;
        }


        for (int i = 0; i < atoms.length; i++) {
            if (lables[i] == 0) {
                buffer.append('-');
            } else
                buffer.append(lables[i]);
        }


        return buffer.toString();

    }

    public void debug() throws Exception {

        Finder finder = new RaphaelFinder(features);

        if (finder.isTRs())
            isTRs = true;

        contentBuilder.append(this.getOutConsole(finder, "RAPHEAL"));

    }


    public void findRepeat() {
        // Phase I. Short TRs detection.
//        ShortTRsDetector detector = new ShortTRsDetector();
//        List<LocationMatch> listFromDetector = detector.getRepeatLocations(features.getStrSeqAp());
//
//        if (listFromDetector.size() > 0)
//            combineScore.setCaScore(1);

        // Phase II. in technical report.


        /**
         * SETUP finders
         */


        List<Finder> finders = new ArrayList<Finder>();
        // conformational alphabets
        ShortTRsFinder shortTRsFinder = new ShortTRsFinder(features);
        finders.add(shortTRsFinder);

        // conformational alphabets finders
        Finder trustfinder = new TrustFinder(features);
        Finder treksfinder = new TReksFinder(features);
        Finder treksSADBFinder = new TReksSADBFinder(features);
        Finder treksOp1Finder = new TReksOp1Finder(features);
        Finder treksOp2Finder = new TReksOp2Finder(features);
        Finder trustCombinedFinder = new TrustCombinedFinder(features);

        finders.add(trustfinder);
        finders.add(treksfinder);
        finders.add(treksSADBFinder);
        finders.add(treksOp1Finder);
        finders.add(treksOp2Finder);
        finders.add(trustCombinedFinder);

        // vector method
        Finder vectorfinder = new VectorFinder(features);
//        Finder vector2ndfinder = new Vector2ndFinder(features);
        Finder vetorsSignalFinder = new VectorSignalFinder(features);

        finders.add(vectorfinder);
//        finders.add(vector2ndfinder);
        finders.add(vetorsSignalFinder);

        // contact map method
        Finder cmFinder = new CMFinder(features);
        Finder pdpFinder = new PDPFinder(features);
        Finder rmsdSignalFinder = new RMSDSignalFinder(features);
        finders.add(cmFinder);
        finders.add(pdpFinder);
        finders.add(rmsdSignalFinder);


        /**
         * setup raphael base finder
         */
        List<Region> regions = this.buildRaphaelFinder(180);
        if (atoms.length < 300) {
            regions.addAll(this.buildRaphaelFinder(atoms.length - 1));
            regions.addAll(this.buildRaphaelFinder((atoms.length - 1) / 2));
            regions.addAll(this.buildRaphaelFinder((atoms.length - 1) / 3));
        }
        for (Region region : regions) {
            RaphaelFinder rapFinder = new RaphaelFinder(features);
            rapFinder.setRegionScan(region);
            finders.add(rapFinder);
        }

        /**
         * end
         */

//        Finder raFinder = new RaphaelFinder(features);
//        finders.add(raFinder);

        // for RMSD finder
//        Finder rmsdFinder = new RMSDFinder(features);
//        finders.add(rmsdFinder);


        /**
         * setup tm score module based on distribution
         * of secondary structure elements
         * from ver 1.1.0
         */
        List<Region> lstRegions = RMSDSubFinder.generateRegions(strSS);
        List<Integer> listPosL = new ArrayList<Integer>();
        for (int i = 2; i < lstRegions.size() / 2; i = i + 1) {
            listPosL.add(i);
        }
        List<List<Integer>> subSets = Lists.partition(listPosL, 1);
        for (List<Integer> subs : subSets) {
            RMSDSubFinder finderR = new RMSDSubFinder(features);
            finderR.setListWinsize(subs);
            finderR.setLstRegions(lstRegions);
            finders.add(finderR);
        }

        /**
         * setup tm score module by using greedy algorithm
         * from ver 1.1.2
         */
        if(atoms.length<300) {
            int maxLenScan = atoms.length / 2;
            int minLenScan = 18;
            if (maxLenScan > 60)
                maxLenScan = 60;
            List<Integer> lstLens = RMSDSubLenFinder.generateLens(minLenScan, maxLenScan);
            List<List<Integer>> subLenSets = Lists.partition(lstLens, 1);
            for (List<Integer> subs : subLenSets) {
                RMSDSubLenFinder finderR = new RMSDSubLenFinder(features);
                finderR.setListWinsize(subs);
                finders.add(finderR);
            }

        }
        // ligand based method
        Finder ligandFinder = new LigandBaseFinder(features);
        finders.add(ligandFinder);

        // CE-Symm method based on internal symmetry of protein structure
        Finder ceSymmFinder = new CESymmFinder(features);
        finders.add(ceSymmFinder);

        /**
         * all finders are run concurrently
         */
        ExecutorService eservice = Executors
                .newFixedThreadPool(nrOfProcessors);
        List<Future> futuresList = new ArrayList<Future>();
        if (atoms.length >= 10 && atoms.length < 1000) {
            for (Finder finder : finders) {
                futuresList.add(eservice.submit(new JobFinder(finder)));
            }
        }
        eservice.shutdown();

        if (DEBUG_PRINT_TRACES)
            System.out.println("End finding");
        // System.out.println("end");
        // end process

        // finders.clear();
        List<Finder> findersOut = new ArrayList<Finder>();
        for (Future future : futuresList) {
            try {
                findersOut.add((Finder) future.get());
                //Finder finder = (Finder) future.get();
            } catch (InterruptedException e) {
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }

        //System.out.println("FINISHED");
        isTRs = false;
        //this.extendAllFinder(finders);
//        System.out.println("combine score...");
        //String strCombinedMesg = this.combineAllFinders(finders, listFromDetector);
        List<Repeat> allTRs = new ArrayList<Repeat>();
        for (Finder finder : findersOut) {
            // get score here
            CombineScore cScore = finder.getCombineScore();
            combineScore.setTmScore(Math.max(cScore.getTmScore(), combineScore.getTmScore()));
            combineScore.setvScore(Math.max(cScore.getvScore(), combineScore.getvScore()));
            combineScore.setlScore(Math.max(cScore.getlScore(), combineScore.getlScore()));
            combineScore.setCeScore(Math.max(cScore.getCeScore(), combineScore.getCeScore()));
            combineScore.setPsimScore(Math.max(cScore.getPsimScore(), combineScore.getPsimScore()));
            combineScore.setSigScore(Math.max(cScore.getSigScore(), combineScore.getSigScore()));
            combineScore.setCaScore(Math.max(cScore.getCaScore(), combineScore.getCaScore()));
            combineScore.setRapScore(Math.max(cScore.getRapScore(), combineScore.getRapScore()));
            if (finder.isTRs()) {
                List<Repeat> repeats = finder.getRepeats();
                try {
                    for (Repeat repeat : repeats) {
                        // convert to pairAlign
                        repeat.setFinderName(finder.getName());
                        if (repeat.getScore() >= ThresholdConfig.QA_THRES && !allTRs.contains(repeat))
                            allTRs.add(repeat);
                    }
                } catch (Exception ex) {
                    logger.error(ex.toString());
                }

            }
        }

        isTRs = false;
        /**
         *  make a decision about repeats using SVM
         */
        //
        // load SVM model
        this.svmScore = SVMWapper.getInstance().getProbEstimatesScore(combineScore.toLibSVMFormat());
        //         classify
        if (svmScore > 0.5) isTRs = true;
        else isTRs = false;


//        if (listFromDetector.size() > 0 || allTRs.size() > 0)
//            isTRs = true;
        //cluster
//        System.out.println("clustering....");
        List<ClusterRepeat> clusters = new ArrayList<ClusterRepeat>();
        if (allTRs.size() > 0) {
            clusters = new ClusterLocation().clusterBaseOnPredictedRegion(allTRs, atoms.length);
            // build predict length base on all TRs
        }

        // display mode
        if (mode.equals(FinderMode.HTML)) {
            // html output
            // contentBuilder.append("<table><tr>");
            contentBuilder.append("<table>");
//            for (Finder finder : finders) {
//                contentBuilder.append(HtmlUtils.getOutHtml(finder,
//                        finder.getName(), atoms));
//            }
//            contentBuilder.append(HtmlUtils.getOutHtmlForTRs(clusters, listFromDetector, atoms));
            contentBuilder.append("</table>");
            this.save();
        } else if (mode.equals(FinderMode.WEBSERVER)) {
            // get output console
            //contentBuilder.append(strCombinedMesg);
            for (Finder finder : finders) {
                contentBuilder.append(this.getOutConsole(finder,
                        finder.getName()));
            }

        } else if (mode.equals(FinderMode.CONSOLE)) {

            // ranking
            try {

                //for short detection.

//                List<Repeat> shortRepeats = detector.convertToRepeats(listFromDetector);
//                for (Repeat repeat : shortRepeats) {
//                    repeat.setFinderName("PhaseI");
//                    repeat.setCluster("clus0");
//                    contentBuilder.append(this.buildOutputOneRepeat(repeat) + "\n");
//                }

                int indexG = 1;
                int noRepeats = 0;
//                for (Repeat top : allTRs) {
//                    top.setCluster("clus" + indexG);
//                    contentBuilder.append(this.buildOutputOneRepeat(top) + "\n");
//                    noRepeats++;
//                }

//                System.out.println("clustering done!");
                for (ClusterRepeat c : clusters) {
                    // System.out.println("-----");
//                    c.sortScoreDESC();
                    //get top
                    //Repeat top = c.getRepeats().get(0);
//                    int limit = 5;
                    int clu_index = 0;
                    for (Repeat top : c.getRepeats()) {
                        if (clu_index == 0) {
                            top.setCluster("cl" + indexG + "_selected");
                        } else
                            top.setCluster("cl" + indexG);
                        contentBuilder.append(this.buildOutputOneRepeat(top) + "\n");
                        this.nrMotifs++;
                        noRepeats++;
                        clu_index++;
//                        limit--;
//                        if (limit < 0)
//                            break;
                    }
                    indexG++;
                }

                if (noRepeats == 0)
                    isTRs = false;

            } catch (Exception ex) {
                logger.error(ex.toString());
            }

        }

        if (DEBUG_PRINT_TRACES)
            System.out.println("DONE ALL STEP");

    }

    private List<Region> buildRaphaelFinder(int L) {

        int shift = L / 2;
        //int L = 180;
        List<Region> regions = new ArrayList<Region>();
        // extract repeat candidates
        for (int i = 0; i < atoms.length - L; i = i + shift) {
            int start = i;
            int end = i + L - 1;
            if (end >= atoms.length - 1)
                end = atoms.length - 1;
            Region region = new Region(start, end);
            if (region.size() > L / 2) {
                regions.add(region);
            }
        }
//        System.out.println(regions);

        return regions;
    }


    public String buildOutputOneRepeat(Repeat top) {
        int start = top.getStart();
        int end = top.getEnd();
        // pdbId_chain Finder avgLeng nUnits RL

        int seqStart = start;
        int seqEnd = end;
        if (!outputRelativeScore) {
            seqStart = atoms[start].getGroup()
                    .getResidueNumber().getSeqNum();
            seqEnd = atoms[end].getGroup()
                    .getResidueNumber().getSeqNum();
        }

        StringBuffer unitsBuffer = new StringBuffer();
        for (RepeatContent unit : top.getRepeats()) {
            int seqS = unit.getStart();
            int seqE = unit.getEnd();
            if (!outputRelativeScore) {
                seqS = atoms[unit.getStart()].getGroup()
                        .getResidueNumber().getSeqNum();
                seqE = atoms[unit.getEnd()].getGroup()
                        .getResidueNumber().getSeqNum();

            }

            unitsBuffer.append(seqS + "-" + seqE + ";");
        }
        String unitsStr = unitsBuffer.toString();
        unitsStr = unitsStr.substring(0, unitsStr.length() - 1);
        String meg = "TRs";
        meg = top.getFinderName();
        return (pdbCode + "_" + pdbChain + "\t" + meg + "\t" + top.getCluster()
                + "\t"
                + top.getRepeats().size() + "\t"
                + NumberFormatUtils.format(top.getAvgLength()) + "\t" + seqStart
                + "-" + seqEnd + "\t" + unitsStr + "\t" + NumberFormatUtils.format(top.getScore()) + "\t" + NumberFormatUtils.format(top.getRankScore()))
                ;
    }

    private void save() {
        // save to output HTML
        String fileDir = "output/" + pdbCode.toLowerCase() + "" + pdbChain
                + ".out";

        DataIO.writeToFile(contentBuilder.toString(), fileDir);

        // save to ouput combined

        String fileDirCombined = "output/" + pdbCode.toLowerCase() + ""
                + pdbChain + ".out.combined";

        DataIO.writeToFile(combinedBuilder.toString(), fileDirCombined);
    }

    public String getOutConsole(Finder finder, String program) {

        List<Repeat> repeats = finder.getRepeats();
        StringBuilder builder = new StringBuilder();
        int type = 1;
        if (repeats.size() > 0 && finder.isTRs()) {
            builder.append("<td valign='top' bgcolor='green' > ");
            builder.append("yes");
        } else {
            builder.append("<td valign='top' bgcolor='lightgray' >");
            builder.append("no");
        }
        builder.append("</td>");

        return builder.toString();

    }

    int nrMotifs = 0;
    double percentCoverage = 0.0;


    public String getCompareResultsBetweenEachMethos(Row row) {

        // tag builder
        String tag = "";
        if (row.getByEyeTot() == 1)// repeat
            tag = "<td valign='top' bgcolor='green' >1</td>";
        else if (row.getByEyeTot() == 2)// not repeat in 3D
            tag = "<td valign='top' bgcolor='red' >2</td>";
        else
            // other
            tag = "<td>" + row.getByEyeTot() + "</td>";
        //

        String conclusion = "";
        if (this.isTRs)// repeat
            conclusion = "<td valign='top' bgcolor='green' >1</td>";
        else
            // not repeat in 3D
            conclusion = "<td valign='top' bgcolor='red' >2</td>";

        String annotation = "<td>" + new SADBAnalysis().annotate(this.strSeqAp)
                + "</td>";

        return "<tr><td>" + pdbCode.toLowerCase() + "_" + pdbChain + "</td>"
                + tag + conclusion + annotation + contentBuilder.toString()
                + "</tr>";
    }

    public String getConsoleOutput() {

        String annotation = new SADBAnalysis().annotate(this.strSeqAp);
        if (annotation.equals(""))
            annotation = "none";

        return pdbCode.toLowerCase() + "_" + pdbChain + "\t" + this.nrMotifs
                + "\t" + annotation + "\t" + atoms.length;
    }

    public String getConsoleOutputDetails() {
        return contentBuilder.toString();
    }

    public boolean isRepeat() {
        return this.isTRs;
    }

    public int getNrMotifs() {
        return this.nrMotifs;
    }

    public double getPercentCoverage() {
        return this.percentCoverage;
    }

    public String getOutput() {
        output = "";
        if (this.isTRs) {
            output = "Chain:"
                    + pdbChain
                    + "  repeat found!  <a class='viewIcon' href='result.php?code="
                    + this.pdbCode + "&chain=" + this.pdbChain + "'></a>";
        }

        return output;
    }


    public double getSvmScore() {
        return svmScore;
    }

    public String getTAPOScore() {
        return combineScore.toString() + "|SVM-Score=" + NumberFormatUtils.format(svmScore);
    }

    public String getStrSeq() {

        return this.strSeq;
    }

    public String getStrSeqAp() {
        return this.strSeqAp;
    }

    public String getStrAcc() {
        return this.strAcc;
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

    public String getStrSS() {
        return strSS;
    }

    public void setStrSS(String strSS) {
        this.strSS = strSS;
    }

    public Atom[] getAtoms() {
        return atoms;
    }

    public void setAtoms(Atom[] atoms) {
        this.atoms = atoms;
    }

    public double[] getCmHistos() {
        return cmHistos;
    }

    public void setCmHistos(double[] cmHistos) {
        this.cmHistos = cmHistos;
    }

    public Features getFeatures() {
        return this.features;
    }

    public String getMode() {
        return mode;
    }

    public void setMode(String mode) {
        this.mode = mode;
    }

    public void setNrOfProcessors(int nrOfProcessors) {
        this.nrOfProcessors = nrOfProcessors;
    }

    public boolean isOutputRelativeScore() {
        return outputRelativeScore;
    }

    public void setOutputRelativeScore(boolean outputRelativeScore) {
        this.outputRelativeScore = outputRelativeScore;
    }

}
