package org.cnrs.crbm.lib.multalign;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.ce.CeParameters;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcMain;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcParameters;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.jama.Matrix;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.repeats.Fragement;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.utils.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.xml.sax.SAXException;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MutilAlign {
    private static Logger logger = LoggerFactory.getLogger(MutilAlign.class);
    String pdbCode = "";
    String pdbChain = "";
    String strAligns = "";
    Atom[] atoms = null;

    private static final String GREEN = "FIVYWML";
    private static final String BLUE = "KR";
    private static final String RED = "EDA";
    private static final String BROWN = "P";
    private static final String ORGANGE = "G";
    private static final String YELLOW = "C";

    StringBuilder builder = new StringBuilder();

    final static boolean debug = true;
    Random rand = new Random();
    ReadFasta fasta = new ReadFasta();
    private double qaScore = 0.0;

    private static final String MSA_APP = MsaAPP.THREEDCOM;

    public double getScore1() {
        return score1;
    }

    private double score1 = 0.0;

    public double getScore2() {
        return score2;
    }

    private double score2 = 0.0;

    public double getScore3() {
        return score3;
    }

    public void setScore3(double score3) {
        this.score3 = score3;
    }

    private double score3 = 0.0;

    public double getScore4() {
        return score4;
    }

    public void setScore4(double score4) {
        this.score4 = score4;
    }

    private double score4 = 0.0;


    public double getQAScore() {
        return this.qaScore;
    }

    public MutilAlign(String code, String chain, String align) {

        this.pdbCode = code;
        this.pdbChain = chain;
        this.strAligns = align;

        try {
            // this.buildAlign();
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            this.mutilAligns(repeatFinder, false);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public MutilAlign(String code, String chain, String align, boolean isConsoleMode) {

        this.pdbCode = code;
        this.pdbChain = chain;
        this.strAligns = align;

        try {
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            this.mutilAligns(repeatFinder, isConsoleMode);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public MutilAlign(String fileDir, String code, String chain, String align) {
        this.pdbCode = code;
        this.pdbChain = chain;
        this.strAligns = align;
        try {
            // this.buildAlign();
            RepeatFinder repeatFinder = new RepeatFinder(fileDir, pdbCode, pdbChain, 0);
            this.mutilAligns(repeatFinder, false);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public MutilAlign(Features features, String align) {
        this.strAligns = align;
        try {
            // this.buildAlign();
            this.mutilAlignsWithScore(features);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public MutilAlign() {

    }

    public static void main(String[] args) {

        // this come from absolute index
        // MutilAlign mutilalign = new MutilAlign("1fbl","A","285-329;330-377;378-426;427-464");

        // this is relative index
//        MutilAlign mutilalign = new MutilAlign("1jah", "A",
//                "48-73;76-103;107-136;138-163");
//        MutilAlign mutilalign = new MutilAlign("1g5q", "A",
//                "65-80;81-98;99-113;123-143");

//        MutilAlign mutilalign = new MutilAlign("2agm", "A",
//                "2-18;19-36;37-55");

//        MutilAlign mutilalign = new MutilAlign("1geq", "B",
//                "4-33;34-81;82-100;110-133;134-157;158-192;193-213");

//        MutilAlign mutilalign = new MutilAlign("1a4y", "A",
//                "30-57;58-86;87-114;115-143;144-171;172-200;201-228;229-257;258-285;286-314;315-342;343-371;372-399;400-428;429-454");

//        MutilAlign mutilalign = new MutilAlign("input/rosmann_model.pdb", "pdbID", " ",
//                "0-25;27-52;62-87;88-113;114-139;153-178;179-204");
//        MutilAlign mutilalign = new MutilAlign("1cgd", "A",
//                "0-2;3-5;6-8;9-11;12-14;15-17;18-20;21-23;24-26");


//        MutilAlign mutilalign = new MutilAlign("1efk", "A",
//                "355-380;381-406;407-430");
//        MutilAlign mutilalign = new MutilAlign("2bib", "A",
//                "203-229;230-254;255-287;288-311;312-334;335-355;356-375;376-395;396-415;416-435;436-455;456-475;476-492;493-512");
        MutilAlign mutilalign = new MutilAlign("2bib", "A", "312-355;356-394;395-431;432-451;452-471;472-491;492-512;513-530");
        DataIO.writeToFile(mutilalign.getOutputHtml(),
                "C:/Users/pdoviet/Desktop/test.html");
        // System.out.println(mutilalign.getOutputHtml());

    }

    /**
     * This class aim to build the contact pattern of contact map of each repeat units with other parts of
     * one protein structure. In one repeat unit region, the position of one residue in contact with other residues outside
     * of this repeat region will bet set to 1. Otherwise, there is no contact with other part of structure, the
     * position of pattern will be set to 0.
     * -------[--------]--[---------]----------
     * |_______|
     * in contact
     * 0000000000100100000001111001010000000000
     *
     * @param cmHistos
     * @param cmapsList
     * @param startRegion
     * @param endRegion
     * @param units
     * @return
     */
    private String buildCMPattern(double[] cmHistos, Map<Integer, List<Integer>> cmapsList, int startRegion, int endRegion, String[] units) {


        StringBuffer bufferCM = new StringBuffer();
        int L = endRegion - startRegion + 1;
        int[] pattern = new int[cmHistos.length];
        for (String unit : units) {
            int start = Integer.parseInt(unit.split("-")[0]);
            int end = Integer.parseInt(unit.split("-")[1]);
            for (int i = start; i <= end; i++) {
                int isContactOutSide = 0;
                if (cmapsList.containsKey(i)) {
                    List<Integer> listCM = cmapsList.get(i);
                    for (int j : listCM) {
                        if (j <= start || end <= j) {
                            isContactOutSide = 1;
                            break;
                        }
                    }
                }
                pattern[i] = isContactOutSide;
            }
        }

        for (int i = 0; i < cmHistos.length; i++) {
            bufferCM.append(pattern[i]);
        }
//        for (int i = 0; i < pattern.length; i++) {
//            bufferCM.append(pattern[i]);
//        }
//        for (int i = endRegion + 1; i < cmHistos.length; i++) {
//            if (i <= cmHistos.length - 1)
//                bufferCM.append(0);
//        }
//        System.out.println(bufferCM.toString());
//        if(bufferCM.toString().length()==cmHistos.length){
//            System.out.println("OK");
//        }


//        for (int i = 0; i < cmHistos.length; i++) {
//            if (startRegion <= i && i <= endRegion) {
//
//                if (cmapsList.containsKey(i)) {
//                    List<Integer> listCM = cmapsList.get(i);
//
//                    for (int j : listCM) {
//                        if (startRegion <= j && j <= endRegion) {
//                            bufferCM.append(1);
//                        } else
//                            bufferCM.append(0);
//                    }
//
//                } else
//                    bufferCM.append(0);
//
//            } else {
//                bufferCM.append(0);
//            }
//        }

        return bufferCM.toString();


    }

    private void mutilAlignsWithScore(Features features) throws Exception {

        // String pdbFile = PdbTools.downloadPDB(pdbCode);

        pdbCode = features.getPdbCode().trim();
        pdbChain = features.getPdbChain().trim();

        atoms = features.getAtoms();
        if (atoms.length < 10)
            throw new StructureException("The number of atoms is small!!!");

        String strSeqs = features.getStrSeq();
        String strCF16 = features.getStrSeqAp();
        //String strCF16 = repeatFinder.get
        String[] strFras = strAligns.split(";");
        int startSegment = Integer.parseInt(strFras[0].split("-")[0]);
        int endSegment = Integer.parseInt(strFras[strFras.length - 1].split("-")[1]);
        String strCM = this.buildCMPattern(features.getCmHistos(), features.getCmapsList(), startSegment, endSegment, strFras);


//        strCM = bufferCM.toString();

        //String strCM = repeatFinder.getCmPattern();
        //String strAcc = repeatFinder.getStrAcc();
        //String strSS = repeatFinder.getStrSS();
        // List<Atom[]> frags = new ArrayList<Atom[]>();


        try {

//            bufferERR.append(strSeqs);
//            bufferERR.append(strAcc);
//            bufferERR.append(strSS);
//            bufferERR.append(strCM);

            MsaOuput ouput = this.getMSA(strFras);
            List<List<Integer>> multiples = ouput.getMsa();
            if (multiples.size() != strFras.length) {
                multiples = this.getMSASequenceMethod(strFras, strSeqs).getMsa();
            }

            if (multiples.size() > 0) {
                TSim tsim = new TSim();
                double tsimScore = 0.0;
                //List<String> msaSeqs = this.        //System.out.println(strSeqs);
//                (multiples, strSeqs);

                List<String> msaCF16A = this.convertToMSA(multiples, strCF16);
                // List<String> msaAcc = this.convertToMSA(multiples, strAcc);
                //List<String> msaSS = this.convertToMSA(multiples, strSS);
                List<String> msaCM = this.convertToMSA(multiples, strCM);

                if (debug) {
//                StringBuffer bf = new StringBuffer();
//
//                for (String s : msaCF16A) {
//                    bf.append(s + "\n");
//                }
//                System.out.println(bf.toString());
//                DataIO.writeToFile(bf.toString(), Dir.MUSTANG_TMP_DIR + "/debug_output/" + repeatFinder.getPdbCode() + "" + repeatFinder.getPdbChain() + "_" + startSegment + "_" + endSegment + ".msa");
//
                }
                //tsimScore += tsim.getTsimScore(msaSeqs);

//            System.out.println(tsim.getTsimScore(msaSeqs));
//            System.out.println(tsim.getTsimScore(msaCF16A));
                //tsimScore += tsim.getTsimScore(msaAcc);
                //tsimScore += tsim.getTsimScore(msaSS);
                this.score3 = tsim.getTsimScore(msaCF16A);
                this.score2 = tsim.getTsimScore(msaCM);

                tsimScore = (this.score2 + this.score3) / 2;

//                if (MSA_APP.equals(MsaAPP.MUSTANG))
//                else if (MSA_APP.equals(MsaAPP.THREEDCOM))
                double stSIM = this.getSTsimScore(strFras, multiples, atoms);
                //stSIM = ouput.getAvgTmScore();

//            System.out.println(stSIM + ": " + tsimScore);
                this.qaScore = (stSIM + tsimScore) / 2;
                this.score1 = stSIM;
                this.score4 = ouput.getAvgTmScore();

                //System.out.println(QA);
            }

        } catch (Exception ex) {
//             ex.printStackTrace();
            //logger.error(bufferERR.toString());
            logger.error(features.getPdbCode() + features.getPdbChain() + " meg: " + ex.getMessage(), ex);
        }


    }

    private MsaOuput getMSASequenceMethod(String[] units, String proteinSeq) {
        ConcurrencyTools.setThreadPoolSingle();
        MsaOuput ouput = new MsaOuput();
        List<List<Integer>> multiples = new ArrayList<List<Integer>>();
        try {
            //convert to sequence
            List<String> seqs = new ArrayList<String>();
            for (String strF : units) {
                int start = Integer.parseInt(strF.split("-")[0]);
                int end = Integer.parseInt(strF.split("-")[1]);
                seqs.add(proteinSeq.substring(start, end + 1));
            }
            List<ProteinSequence> lst = new ArrayList<ProteinSequence>();
            int seqName = 0;
            for (String seq : seqs) {
                ProteinSequence sequence = new ProteinSequence(seq);
                sequence.setAccession(new AccessionID("seq" + seqName));
                lst.add(sequence);
                seqName++;
            }

            Profile<ProteinSequence, AminoAcidCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
//            System.out.printf("Clustalw:%n%s%n", profile);
            ConcurrencyTools.shutdown();
            List<AlignedSequence<ProteinSequence, AminoAcidCompound>> alSeq = profile.getAlignedSequences();
            Map<String, String> mapSeq = new HashMap<String, String>();
            for (Sequence<AminoAcidCompound> seq : alSeq) {
                mapSeq.put(seq.getAccession().getID(), seq.getSequenceAsString());
            }
//            String[] msa = profile.toString().split("\n");
            //System.out.println(msa.split("\n")[0]);
//            System.out.println("--" + msa[0] + "--");
//            int msaLeng = msa[0].length();

            int msaLeng = alSeq.get(0).getSequenceAsString().length();
            //System.out.println(msaLeng);
            if (msaLeng > 0) {
                for (int index = 0; index < units.length; index++) {
                    int startRegion = Integer.parseInt(units[index].split("-")[0]);
                    int endRegion = Integer.parseInt(units[index].split("-")[1]);
                    List<Integer> list = new ArrayList<Integer>();
                    String unit = mapSeq.get("seq" + index); // msa[index];
                    int k = 0;
                    for (int j = 0; j < msaLeng; j++) {
                        char AA = unit.charAt(j);
                        int pos = -1;
                        if (AA != '-') {
                            pos = startRegion + k;
//                            if(pos>endRegion)
//                                pos = endRegion;
//                            if (pos == 125) {
//                                //System.out.println(lst);
//                                System.out.println("------------");
//                                System.out.println(seqs);
//                                System.out.println(profile.toString());
//                                System.out.println(units[index] + ":" + index);
//                                System.out.println(unit);
//                                System.out.println(start);
//                            }
                            k++;
                        } else {
                            pos = -1;
                        }
                        list.add(pos);
                    }
                    multiples.add(list);
                }
            }

        } catch (Exception ex) {
            //do nothing :)
//            ex.printStackTrace();
//            System.out.println(proteinSeq);
            logger.warn("MutilAlign", ex);
        }
        ouput.setMsa(multiples);
        return ouput;
    }


    /**
     * BUG>......
     *
     * @param units
     * @return
     */
    private MsaOuput getMSA_NextGen(String[] units) {
        MsaOuput ouput = new MsaOuput();
        List<List<Integer>> multiples = new ArrayList<List<Integer>>();
        try {
            List<Atom[]> atomArrays = new ArrayList<Atom[]>();
            for (String strF : units) {
                int start = Integer.parseInt(strF.split("-")[0]);
                int end = Integer.parseInt(strF.split("-")[1]);
                if (end >= atoms.length)
                    end = atoms.length - 1;
                atomArrays.add(Fragement.getFragementsofAtomsWithClone(atoms, start, end));

            }
//            ConcurrencyTools.setThreadPoolSingle();
            //Here the multiple structural alignment algorithm comes in place to generate the alignment object

//            List<String> lstXXX = new ArrayList<String>();
//            for (String u : units) {
//                lstXXX.add(u);
//            }
//            System.out.println(lstXXX);
            MultipleMcMain algorithm = new MultipleMcMain(new CeMain());
            MultipleMcParameters params = (MultipleMcParameters) algorithm.getParameters();
//            params.setMinBlockLen(10);
//            params.setGapExtension(20.0);
            params.setNrThreads(1);
            MultipleAlignment result = algorithm.align(atomArrays);
            result.getEnsemble().setStructureNames(Arrays.asList(units));
//            ConcurrencyTools.shutdown();
//            System.out.println("done!");
            //System.out.println(result.getScore("AvgTM-score"));
//            System.out.println(MultipleAlignmentWriter.toFASTA(result));
            List<List<Integer>> multiples_tmp = MSAWriter.tapoMsaFormat(result);
            // convert to
            for (int i = 0; i < multiples_tmp.size(); i++) {
                int start = Integer.parseInt(units[i].split("-")[0]);
                List<Integer> lst = new ArrayList<Integer>();
                for (int val : multiples_tmp.get(i)) {
                    lst.add(val + start);
                }
                multiples.add(lst);

            }

        } catch (Exception ex) {
//            ex.printStackTrace();
            logger.error(ex.toString());
//            List<String> lst = new ArrayList<String>();
//            for (String u : units) {
//                lst.add(u);
//            }
//            System.err.println(lst);
//            ex.printStackTrace();
//            System.exit(1);
        }
        ouput.setMsa(multiples);
        return ouput;
    }


    private MsaOuput getMSA(String[] units) {
        MsaOuput ouput = new MsaOuput();
        List<List<Integer>> multiples = new ArrayList<List<Integer>>();
        int startSegment = Integer.parseInt(units[0].split("-")[0]);
        int endSegment = Integer.parseInt(units[units.length - 1].split("-")[1]);
        String fileExt = "msa_tmp_" + pdbCode + pdbChain.trim() + "_" + startSegment + "_" + endSegment + "_"
                + MD5Utils.MD5(rand.nextInt(10000000) + "");
//        String fileExt = "msa_tmp_" + pdbCode + pdbChain + "_" + startSegment + "_" + endSegment ;

        String tmpDir = Dir.DCOMB_TMP_DIR + "/" + fileExt;
        File theTmpDir = new File(tmpDir);
        if (!theTmpDir.exists()) {
            // make a directory before save file
            theTmpDir.mkdir();
        }
        int i = 0;
        try {
            StringBuffer buffer = new StringBuffer();
            if (MSA_APP.equals(MsaAPP.MUSTANG))
                buffer.append("> " + tmpDir + "\n");
            for (String strF : units) {
                int start = Integer.parseInt(strF.split("-")[0]);
                int end = Integer.parseInt(strF.split("-")[1]);
                if (end >= atoms.length)
                    end = atoms.length - 1;
                PdbTools.writePDB(atoms, pdbCode, pdbChain, start, end, theTmpDir
                        + "/" + i + ".pdb");
                if (MSA_APP.equals(MsaAPP.MUSTANG))
                    buffer.append("+" + i + ".pdb\n");
                else
                    buffer.append(tmpDir + "/" + i + ".pdb\n");
                i++;
            }

            DataIO.writeToFile(buffer.toString(), tmpDir + "/pdbs.in");
            String outputFile = "";
            if (MSA_APP.equals(MsaAPP.MUSTANG))
                outputFile = tmpDir + "/msa";
            else
                outputFile = tmpDir + "/pdbs.in.ali";

            /**
             *  call MSA program
             */
            ouput = calBashshell(tmpDir + "/pdbs.in", outputFile);

            boolean msaFileExis = false;
            Map<String, String> a = null;
            if (MSA_APP.equals(MsaAPP.MUSTANG)) {
                if (new File(outputFile + ".afasta").exists()) {
                    a = fasta.getFastaFileInSpecialCase(outputFile + ".afasta");
                    msaFileExis = true;
                }
                // else throw new Exception(outputFile + ".afasta does not exist!!");
            } else { // 3DCOMB application
                if (new File(outputFile).exists()) {
                    a = fasta.getFastaFileInSpecialCase(outputFile);
                    msaFileExis = true;
                }
                // else throw new Exception(outputFile + " does not exist!!");
            }
//            Entry<String, ProteinSequence> mapEntry = a.entrySet().iterator().next();
//            int msaLeng = mapEntry.getValue().getSequenceAsString().length();

            if (msaFileExis) {
                int msaLeng = a.entrySet().iterator().next().getValue().length();
                if (msaLeng > 0) {
                    for (int index = 0; index < units.length; index++) {
                        String name = tmpDir + "/" + index + ".pdb";
                        if (MSA_APP.equals(MsaAPP.MUSTANG))
                            name = index + ".pdb";
                        //System.out.println(seqs);
                        int start = Integer.parseInt(units[index].split("-")[0]);
                        List<Integer> list = new ArrayList<Integer>();
                        if (a.containsKey(name)) {
                            //String seqs = a.get(name).getSequenceAsString();
                            String seqs = a.get(name.trim());
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
                        } else {
                            for (int j = 0; j < msaLeng; j++) {
                                list.add(-1);
                            }

                        }
                        multiples.add(list);
                    }
                }
            }

        } catch (Exception ex) {
            //ex.printStackTrace();
            logger.error(ex.getMessage(), ex);
        }
        try {
            FileUtils.deleteDirectory(theTmpDir);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            // e.printStackTrace();
        }

        ouput.setMsa(multiples);
        return ouput;
    }


    private void mutilAligns(RepeatFinder repeatFinder, boolean isConsoleMode) throws Exception {
        // String pdbFile = PdbTools.downloadPDB(pdbCode);
        //RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        // Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        atoms = repeatFinder.getAtoms();
        if (atoms.length < 10)
            throw new StructureException("The number of atoms is small!!!");
        String strAcc = repeatFinder.getStrAcc();
        String strSS = repeatFinder.getStrSS();
        String strSeqs = repeatFinder.getStrSeq();
//        System.out.println(strSeqs);
        String[] strFras = strAligns.split(";");
        int startSegment = Integer.parseInt(strFras[0].split("-")[0]);
        int endSegment = Integer.parseInt(strFras[strFras.length - 1].split("-")[1]);
        //String strCM = this.buildCMPattern(repeatFinder.getCmHistos(), repeatFinder.getCmapsList(), startSegment, endSegment);
        String strCF16 = repeatFinder.getStrSeqAp();

        try {

            List<List<Integer>> multiples = new ArrayList<List<Integer>>();
            if (strFras.length > 20 && strFras[0].length() < 10) {
                multiples = this.getMSASequenceMethod(strFras, strSeqs).getMsa();
            } else {
                MsaOuput ouput = this.getMSA(strFras);
                multiples = ouput.getMsa();
            }
            if (multiples.size() != strFras.length) {
                multiples = this.getMSASequenceMethod(strFras, strSeqs).getMsa();
            }
            //List<String> msaCM = this.convertToMSA(multiples, strCM);
//            List<String> msaCF16A = this.convertToMSA(multiples, strCF16);
//            if (debug) {
//
//                StringBuffer bf = new StringBuffer();
//
//                for (String s : msaCF16A) {
//                    bf.append(s + "\n");
//                }
//                System.out.println(bf.toString());
//                //DataIO.writeToFile(bf.toString(), Dir.MUSTANG_TMP_DIR + "/msa_output/" + repeatFinder.getPdbCode() + "" + repeatFinder.getPdbChain() + ".msa");
//            }
//            System.err.println(this.getSTsimScore(strFras, multiples, atoms));
//            this.getSTsimScore(strFras, multiples, atoms);
            if (isConsoleMode) {

                //List<String> msaCM = this.convertToMSA(multiples, strCM);
                List<String> msaSeq = this.convertToMSA(multiples, strSeqs);
                //StringBuffer bf = new StringBuffer();
                for (String s : msaSeq) {
                    builder.append(s + "\n");
                }
//                System.out.println(builder.toString());

            } else

                builderHTML(multiples, strAcc, strSS, strFras, strSeqs);

        } catch (Exception ex) {

            //ex.printStackTrace();
            logger.error(ex.getMessage(), ex);
        }


    }


    public double getSTsimScore(String[] units, List<List<Integer>> msa, Atom[] atoms) {
        double m = units.length;
        int N = 0;
        double score = 0.0;
        double maxTMScore = 0.0;
        try {
            for (int i = 0; i < m - 1; i++) {

                int starti = Integer.parseInt(units[i].split("-")[0]);
                int endi = Integer.parseInt(units[i].split("-")[1]);
                Atom[] atomi = Fragement.getFragementsofAtoms(atoms, starti, endi);

                for (int j = i + 1; j < m; j++) {
                    int startj = Integer.parseInt(units[j].split("-")[0]);
                    int endj = Integer.parseInt(units[j].split("-")[1]);
                    Atom[] atomj = Fragement.getFragementsofAtoms(atoms, startj, endj);

//                    AFPChain afpChain = this.pairAlign(atomi, atomj);
//                    System.out.println(afpChain.getTMScore());
                    try {
                        List<Integer> struct1 = msa.get(i);
                        List<Integer> struct2 = msa.get(j);
                        List<Atom> latom1 = new ArrayList<Atom>();
                        List<Atom> latom2 = new ArrayList<Atom>();
                        for (int index = 0; index < struct1.size(); index++) {
                            if (struct1.get(index) != -1 && struct2.get(index) != -1) {
                                latom1.add((Atom) atoms[struct1.get(index)].clone());
                                latom2.add((Atom) atoms[struct2.get(index)].clone());
                            }
                        }
                        Atom[] ca1aligned = latom1.toArray(new Atom[latom1.size()]);
                        Atom[] ca2aligned = latom2.toArray(new Atom[latom2.size()]);

                        //Superimpose
                        if (ca1aligned.length > 0) {
                            SVDSuperimposer svd =
                                    new SVDSuperimposer(ca1aligned, ca2aligned);
                            Matrix matrix = svd.getRotation();
                            Atom shift = svd.getTranslation();
                            for (Atom a : ca2aligned) {
                                Calc.rotate(a, matrix);
                                Calc.shift(a, shift);
                            }

                            double tm = AFPChainScorerAlter.getTMScore(ca1aligned, ca2aligned, atomi.length, atomj.length);

//                            System.out.println(tm);
                            maxTMScore = Math.max(maxTMScore, tm);
                            score += tm;
                        }
                    } catch (Exception e) {
                        // e.printStackTrace();
                    }

                    N++;
                }

            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        //if (N == 0) return 0.0;
        //else return score / N;
        return maxTMScore;

    }


    class MsaOuput {
        private double avgTmScore = 0.0;

        public double getnCore() {
            return nCore;
        }

        public void setnCore(double nCore) {
            this.nCore = nCore;
        }

        public double getAvgTmScore() {
            return avgTmScore;
        }

        public void setAvgTmScore(double avgTmScore) {
            this.avgTmScore = avgTmScore;
        }

        private double nCore = 0;

        public double getRmsd() {
            return rmsd;
        }

        public void setRmsd(double rmsd) {
            this.rmsd = rmsd;
        }

        private double rmsd;

        public List<List<Integer>> getMsa() {
            return msa;
        }

        public void setMsa(List<List<Integer>> msa) {
            this.msa = msa;
        }

        List<List<Integer>> msa = new ArrayList<List<Integer>>();

    }

    private MsaOuput calBashshell(String inputFile, String outputFile) {

        // new version from 2.0.4.2
        MsaOuput ouput = new MsaOuput();

        String shellCommand = Dir.DCOMB_EXECUTABLE + " " + inputFile;
        if (MSA_APP.equals(MsaAPP.MUSTANG))
            shellCommand = Dir.MUSTANG_EXECUTABLE + " -f " + inputFile
                    + " -o " + outputFile + " -F fasta -s OFF";

        String strCmdOut = ExecUtils.execShellCmdLinux(shellCommand);
        if (MSA_APP.equals(MsaAPP.THREEDCOM)) {
            String[] lines = strCmdOut.split("\n");
            List<Double> numbers = new ArrayList<Double>();
            for (String line : lines) {
                //System.out.println(line);

                if (line.startsWith("CORE_LEN")) {
                    Pattern p = Pattern.compile("=([^,]*)");
                    Matcher m = p.matcher(line);

                    try {

                        while (m.find()) {
                            //System.out.println(m.group().replaceAll("=", "").trim());
                            numbers.add(Double.parseDouble(m.group().replaceAll("=", "").trim()));
                        }

                    } catch (Exception ex) {
                        logger.error(ex.getMessage(), ex);
                    }
                }
            }


            if (numbers.size() == 3) {
                ouput.setnCore(numbers.get(0));
                ouput.setRmsd(numbers.get(1));
                ouput.setAvgTmScore(numbers.get(2));
            }
        }
        return ouput;

    }

//    private void buildAlign() throws StructureException {
//
//        // String pdbFile = PdbTools.downloadPDB(pdbCode);
//
//        //RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
//
//        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
//        atoms = StructureTools
//                .getAtomCAArray(structure.getChainByPDB(pdbChain));
//        if (atoms.length < 10)
//            throw new StructureException("The number of atoms is small!!!");
//        DSSP dssp = new DSSP(pdbCode);
//
//        String strAcc = dssp.getRelativeACC(pdbChain, atoms);
//        String strSS = dssp.filterSS(dssp.getCombinedSS(atoms, this.pdbChain));
//        String strSeqs = dssp.getAminoAcidSequence(pdbChain);
//        // List<Atom[]> frags = new ArrayList<Atom[]>();
//
//        String[] strFras = strAligns.split(";");
//
//        // if (debug) {
//        // for (String strF : strFras) {
//        // int start = Integer.parseInt(strF.split("-")[0]);
//        // int end = Integer.parseInt(strF.split("-")[1]);
//        // Atom[] frag = Fragement.getFragementsofAtoms(atoms, start, end);
//        // // frags.add(frag);
//        // PdbTools.writePDB(atoms, pdbCode, pdbChain, start, end);
//        //
//        // }
//        // }
//
//        // take a seed
//
//        String seedstrF = strFras[0];
//        int seedstart = Integer.parseInt(seedstrF.split("-")[0]);
//        int seedend = Integer.parseInt(seedstrF.split("-")[1]);
//        Atom[] seed = Fragement.getFragementsofAtoms(atoms, seedstart, seedend);
//
//        int nRows = strFras.length;
//        int nCols = seedend - seedstart + 1;
//        int[][] multi = new int[nRows][nCols];
//        int row = 0;
//        for (int i = seedstart; i <= seedend; i++) {
//            multi[row][i - seedstart] = i;
//        }
//
//        row = 1;
//        for (int i = 1; i < strFras.length; i++) {
//            int start = Integer.parseInt(strFras[i].split("-")[0]);
//            int end = Integer.parseInt(strFras[i].split("-")[1]);
//            Atom[] atomcompare = Fragement.getFragementsofAtoms(atoms, start,
//                    end);
//            AFPChain afpChain = this.pairAlign(seed, atomcompare);
//            Map<Integer, Integer> alignPairs = this.toAlignedPairs(afpChain,
//                    seedstart, start);
//            // int startCore = 1000;
//            // int endCore = 0;
//            for (int j = seedstart; j <= seedend; j++) {
//                if (alignPairs.containsKey(j)) {
//                    multi[row][j - seedstart] = alignPairs.get(j);
//                    // startCore = Math.min(startCore, alignPairs.get(j));
//                    // endCore = Math.max(endCore, alignPairs.get(j));
//                } else
//                    multi[row][j - seedstart] = -1;
//            }
//
//            row++;
//        }
//
//        List<List<Integer>> multiples = AlignUtils.getMultipleAligns(multi,
//                nRows, nCols);
//
//        //builderHTML(multiples, strAcc, strSS, strFras, strSeqs);
//
//    }


    private int getResSeq(int pos) {

        return atoms[pos].getGroup().getResidueNumber().getSeqNum();
    }


    public List<String> convertToMSA(List<List<Integer>> multiples, String strSeqs) throws Exception {
//        System.out.println(strSeqs);
        List<String> msa = new ArrayList<String>();
        for (int k = 0; k < multiples.size(); k++) {
            StringBuffer buffer = new StringBuffer();
            List<Integer> align = multiples.get(k);
            for (Integer index : align) {
                if (index != -1) {
//                    try {
                    buffer.append(strSeqs.charAt(index));

//                    } catch (Exception ex) {
//                        String meg = strSeqs + "\n" + strSeqs.length() + ":" + index;
//                        throw new Exception(meg);
//                    }

                } else {
                    buffer.append("-");
                }
            }
            msa.add(buffer.toString());
        }


        return msa;
    }

    private void builderHTML(List<List<Integer>> multiples, String strAcc,
                             String strSS, String[] strFags, String strSeqs) {
        // get string aligns
        try {
            builder.append("<table border='0' cellspacing='0' cellpadding='2'>\n");
            builder.append("<tbody>\n");

            for (int k = 0; k < multiples.size(); k++) {
                builder.append("<tr>");
                int start = Integer.parseInt(strFags[k].split("-")[0]);
                int end = Integer.parseInt(strFags[k].split("-")[1]);
                List<Integer> align = multiples.get(k);
                int startCore = 1000;
                int endCore = 0;

                for (Integer index : align) {
                    if (index != -1) {
                        startCore = Math.min(startCore, index);
                        endCore = Math.max(endCore, index);
                    }

                }

                // int abudant = k % 5;

                int size = ColorUtils.getInstance().getColors().size();
                int abudant = k % size;

                int red = ColorUtils.getInstance().getColors().get(abudant)
                        .getRed();
                int green = ColorUtils.getInstance().getColors().get(abudant)
                        .getGreen();
                int blue = ColorUtils.getInstance().getColors().get(abudant)
                        .getBlue();
                // String bg_repeat = "bgcolor='rgb(" + red + "," + green + "," +
                // blue
                // + ")'";
                String bg_repeat = "style='width:10px ; background-color:rgb("
                        + red + "," + green + "," + blue + ")'";
                // if (abudant == 0)
                // bg_repeat = "bgcolor='yellow'";
                // else if (abudant == 1)
                // bg_repeat = "bgcolor='green'";
                // else if (abudant == 2)
                // bg_repeat = "bgcolor='red'";
                // else if (abudant == 3)
                // bg_repeat = "bgcolor='blue'";
                // else if (abudant == 4)
                // bg_repeat = "bgcolor='orange'";

                builder.append("<td  valign='center' " + bg_repeat + "></td>");

                builder.append(String
                        .format("<td  valign='center'><font  size='1' face='Courier New,Courier'><font >%s</font></td>",
                                this.getResSeq(start)));


                if (startCore > 0 && startCore >= start && startCore < 1000) {
                    String left = strSeqs.substring(start, startCore);
                    builder.append(String
                            .format("<td pairAlign='right' valign='center'><font face='Courier New,Courier'><font >%s</font></td>",
                                    left.toLowerCase()));
                } else
                    builder.append("<td></td>");
                if (startCore < 1000 && endCore > 0) {
                    for (Integer index : align) {
                        // fix bug;
                        if (index > end || index < start)
                            index = -1;
                        if (index == -1)
                            builder.append(String.format("<td>%s</td>", "-"));
                        else {
                            char AA = strSeqs.charAt(index);
                            //char AA = this.getOneLetter(atoms[index].getGroup());
                            char ACC = strAcc.charAt(index);
                            char struct = strSS.charAt(index);
                            String residuce = "";
                            String bgcolor = "";
                            String fontcolor = "color='" + this.getFontColor(AA) + "'";
                            String underline = "";
                            residuce = "" + AA + "";
                            if (ACC == '0') {
                                bgcolor = "bgcolor='lightgray'";
                                residuce = "<b>" + AA + "</b>";
                            }

                            if (struct == 'H') {
                                underline = "style='border-bottom:2px solid black'";
                            } else if (struct == 'B') {
                                underline = "style='border-bottom:2px solid red'";
                            } else
                                underline = "style='border-bottom:2px solid white'";
                            builder.append(String
                                    .format("<td pairAlign='center' valign='center' "
                                                    + bgcolor
                                                    + " "
                                                    + underline
                                                    + " ><font face='Courier New,Courier'><font "
                                                    + fontcolor + ">%s</font></font></td>",
                                            residuce));

                        }

                    }


                } else {
                    // there are no multiple alignment
//                    String noalign = strSeqs.substring(start, end);
//                    builder.append("<td></td>");

                    for (Integer index : align) {
                        // fix bug;
                        if (index > end || index < start)
                            index = -1;
                        if (index == -1)
                            builder.append(String.format("<td>%s</td>", "-"));
                    }

                }

                if (endCore < 1000 && endCore <= end && endCore > 0) {
                    String right = strSeqs.substring(endCore, end);
                    builder.append(String
                            .format("<td pairAlign='right' valign='center'><font face='Courier New,Courier'><font >%s</font></td>",
                                    right.toLowerCase()));
                } else {
                    builder.append("<td></td>");
                }
                builder.append(String
                        .format("<td  valign='center'><font  size='1' face='Courier New,Courier'><font >%s</font></td>",
                                this.getResSeq(end)));

                builder.append("</tr>\n");

            }

            // System.out.println(linkView);
            builder.append("</tbody>\n");
            builder.append("</table>\n");
            builder.append("\n");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private char getOneLetter(Group g) {

        try {
//            Character c = StructureTools.get1LetterCode(g.getPDBName());
//            return c;
            return g.getChemComp().getOne_letter_code().charAt(0);
        } catch (Exception e) {
            return 'X';
        }
    }

    public double alignContactMap(AFPChain afpChain, int startSet1,
                                  int startSet2, double[] cmHisto) {

        double rmsd = 10.0;
        // afpChain.getCa1Length();

        try {
            Map<Integer, Integer> alignPairs = this.toAlignedPairs(afpChain,
                    startSet1, startSet2);

            double square = 0.0;
            int leng = alignPairs.size();
            for (Map.Entry<Integer, Integer> entry : alignPairs.entrySet()) {

                double histo1 = cmHisto[entry.getKey()];
                double histo2 = cmHisto[entry.getValue()];

                if (histo1 > leng)
                    histo1 = 1;
                else
                    histo1 = 0;
                if (histo2 > leng)
                    histo2 = 1;
                else
                    histo2 = 0;

                square += Math.pow(histo1 - histo2, 2);
            }

            rmsd = (double) square / alignPairs.size();

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return rmsd;

    }

    public AFPChain pairAlign(Atom[] atomSet1, Atom[] atomSet2)
            throws StructureException {
        try {
            // configure method
            //ConfigStrucAligParams params = new FatCatParameters();
            //String algorithmName = FatCatFlexible.algorithmName;
            ConfigStrucAligParams params = new CeParameters();
            String algorithmName = CeMain.algorithmName;
            int noRes1 = atomSet1.length;
            int noRes2 = atomSet2.length;
            if (noRes1 > 0) {
//                if ((double) (Math.abs(atomSet2.length - atomSet1.length) / atomSet1.length) > 0.2)
//                    noRes1 = (int) (atomSet2.length + atomSet1.length) / 2;

//                /**
//                 * from version 1.1.2-dev
//                 */
//                if (noRes1 < 20)// change algorithm
//                {
//                    params = new FatCatParameters();
//                    algorithmName = FatCatRigid.algorithmName;
//                }

                StructureAlignment algorithm = StructureAlignmentFactory
                        .getAlgorithm(algorithmName);
                AFPChain afpChain = algorithm.align(atomSet1, atomSet2, params);
                // System.out.println(afpChain.getChainRmsd() + " "
                // + afpChain.getChainLen());
                double tmScore = AFPChainScorerAlter.getTMScore(afpChain, atomSet1,
                        atomSet2);
                if (tmScore <= 0) tmScore = 0;
                afpChain.setTMScore(tmScore);
                return afpChain;
            }
        } catch (Exception e) {
            // System.out.println("set 1: " + atomSet1.length + "  set 2: "
            // + atomSet2.length);
            e.printStackTrace();
        }
        return null;

    }


    public AFPChain pairAlignTMScoreDefault(Atom[] atomSet1, Atom[] atomSet2)
            throws StructureException {
        try {
            // configure method
            //ConfigStrucAligParams params = new FatCatParameters();
            //String algorithmName = FatCatFlexible.algorithmName;
            ConfigStrucAligParams params = new CeParameters();
            String algorithmName = CeMain.algorithmName;
            int noRes1 = atomSet1.length;
            int noRes2 = atomSet2.length;
            if (noRes1 > 0) {


                StructureAlignment algorithm = StructureAlignmentFactory
                        .getAlgorithm(algorithmName);
                AFPChain afpChain = algorithm.align(atomSet1, atomSet2, params);
                // System.out.println(afpChain.getChainRmsd() + " "
                // + afpChain.getChainLen());
                double tmScore = AFPChainScorer.getTMScore(afpChain, atomSet1, atomSet2);
//                double tmScore = AFPChainScorerAlter.getTMScore(afpChain, atomSet1,
//                        atomSet2);

                if (tmScore <= 0) tmScore = 0;
                afpChain.setTMScore(tmScore);
                return afpChain;
            }
        } catch (Exception e) {

            e.printStackTrace();
        }
        return null;

    }

    public Map<Integer, Integer> toAlignedPairs(AFPChain afpChain,
                                                int abPositionA1, int abPositionA2) {
        Map<Integer, Integer> aligns = new HashMap<Integer, Integer>();
        int[][][] optAln = afpChain.getOptAln();
        int[] blockLen = afpChain.getOptLen();
        for (int block = 0; block < afpChain.getBlockNum(); block++) {
            for (int i = 0; i < blockLen[block]; i++) {
                int posA1 = abPositionA1 + optAln[block][0][i];
                int posA2 = abPositionA2 + optAln[block][1][i];
                aligns.put(posA1, posA2);

            }
        }
        return aligns;

    }

    private String getFontColor(char AA) {
        String residuce = AA + "";
        if (GREEN.contains(residuce))
            return "green";
        if (YELLOW.contains(residuce))
            return "yellow";
        if (RED.contains(residuce))
            return "red";
        if (BLUE.contains(residuce))
            return "blue";
        if (BROWN.contains(residuce))
            return "brown";
        if (ORGANGE.contains(residuce))
            return "orange";

        return "black";

    }

    public String getOutput() {

        return builder.toString();
    }

    public String getOutputHtml() {

        return builder.toString();
    }

    public BufferedImage getImage() throws IOException, SAXException {
        // HtmlImageGenerator imageGenerator = new HtmlImageGenerator();
        // imageGenerator.loadHtml(builder.toString());
        // return imageGenerator.getBufferedImage();

        ImageRenderPng imageRender = new ImageRenderPng();
        return imageRender.renderURL(builder.toString());

    }
}
