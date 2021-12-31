package org.cnrs.crbm.lib.dssp;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.secstruc.SecStrucInfo;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.sadb.StrucState;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class DSSP {

    private final static double THRES_ACC = 0.20;
    static Logger logger = LoggerFactory.getLogger(DSSP.class);

    public static void main(String[] args) throws StructureException {

        // String pdbCode = args[0];

        String pdbCode = "2z3s";
        String pdbChain = "A";

        // String pdbFile = PdbTools.downloadPDB(pdbCode);
        // Structure structure = PdbTools.getStructureFromFile(pdbFile);
        //

        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        Atom[] atoms = PdbTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain));
        DSSP dssp = new DSSP(pdbCode, pdbChain);
        System.out.println(dssp.getAminoAcidSequence("A"));

        // DSSP dssp = new DSSP(pdbCode);
        // System.out.println(dssp.getAminoAcidSequence("A"));

//        System.out.println("before:");
//        System.out.println(dssp.getSS("A", atoms));
//        System.out.println("after:");
//        System.out.println(dssp.filterSS(dssp.getSS("A", atoms)));


//        String strAcc = dssp.getRelativeACC(pdbChain, atoms);
//
//        for (int i = 0; i < strAcc.length(); i++) {
//
//            Atom atom = atoms[i];
//
//            if ('0' == strAcc.charAt(i)) {
//                int seqNum = atom.getGroup().getResidueNumber().getSeqNum();
//                // System.out.print(seqNum+"+");
//                System.out.println("color red, /" + pdbCode + "//" + pdbChain + "/" + seqNum);
//            }
//
//
//        }

        //System.out.println(strAcc);
        // DSSP dssp = new DSSP();
        // dssp.parser("output/1LXA.dssp", true);
        // System.out.println(dssp.getAminoAcidSequence("A"));
    }


    public static boolean READ_CLASSIC_FORMAT = false;


    /**
     * Auto insert to dssp database if there is missing this dssp file
     * format of this application will be e.g 1lxaA.dssp
     *
     * @param pdbCode
     * @param pdbChain
     */
    public DSSP(String pdbCode, String pdbChain) {
        this.c = new ArrayList<DSSPDataLine>();

        String dsspFile = Dir.DSSP_LOCAL + "/" + pdbCode.substring(1, 3) + "/"
                + pdbCode + pdbChain + ".dssp";
        // System.out.println(dsspFile);
        boolean isExists = new File(dsspFile).exists();
        if (!isExists) {
            // insert
            DsspLocalDb dsspDb = new DsspLocalDb();
            dsspDb.insertDB(pdbCode, pdbChain);
        }
        this.parser(dsspFile, READ_CLASSIC_FORMAT);

        // mapping here

    }


    /**
     * DSSP with file model.
     *
     * @param fileDir
     * @param pdbCode
     * @param pdbChain
     */
    public DSSP(String fileDir, String pdbCode, String pdbChain) {
        this.c = new ArrayList<DSSPDataLine>();
        DsspLocalDb dsspDb = new DsspLocalDb();
        String dsspFile = Dir.TMP_DIR + "/"
                + pdbCode + pdbChain.trim() + ".dssp";
        dsspDb.generateDSSP(fileDir, dsspFile);
        this.parser(dsspFile, READ_CLASSIC_FORMAT);
        new File(dsspFile).delete();

        // mapping here

    }

    /**
     * DSSP with no processing
     */
    public DSSP() {

    }

    /**
     * Auto insert to dssp database if there is missing this dssp file.
     *
     * @param pdbCode
     */
    public DSSP(String pdbCode) {
        this.c = new ArrayList<DSSPDataLine>();

        String dsspFile = Dir.DSSP_LOCAL + "/" + pdbCode.substring(1, 3) + "/"
                + pdbCode + ".dssp";
        // System.out.println(dsspFile);
        boolean isExists = new File(dsspFile).exists();
        if (!isExists) {
            // insert
            DsspLocalDb dsspDb = new DsspLocalDb();
            dsspDb.insertDB(pdbCode);
        }
        this.parser(dsspFile, READ_CLASSIC_FORMAT);

        // mapping here

    }

    private int a = -1;
    private ArrayList<DSSPDataLine> c = new ArrayList<DSSPDataLine>();

    private byte d = 0;

    public ArrayList<DSSPDataLine> getDSSPDataLines() {
        return this.c;
    }

    public String getAminoAcidSequence(String chainName) {

        StringBuilder builder = new StringBuilder();
        for (DSSPDataLine line : c) {
            if (chainName.equals(line.chainLetter + ""))
                builder.append(line.residueLetter);
        }
        return builder.toString();
    }

    /**
     * Last update 28/03: consider to the contact map - HHHHH <=5 and BB<2 were
     * considered as turn.
     *
     * @author phuongdo
     */
    public String filterSS(String secondaryStructure) {
        StringBuilder builder = new StringBuilder();

        // configuration parameters here
        int N_HELIX = 5;
        int N_STRAND = 1;

        int nres = secondaryStructure.length();
        for (int i = 0; i < nres; i++) {
            char c = secondaryStructure.charAt(i);
            if (c == '-')// turn
                builder.append("-");
            else if (c == 'B') {
                StringBuilder sBuilder = new StringBuilder();
                int l = i;
                char cB = 'B';
                do {

                    sBuilder.append(cB);
                    l++;
                    if (l < nres)
                        cB = secondaryStructure.charAt(l);

                } while (cB == 'B' && l < nres);

                String ss = sBuilder.toString();
                if (ss.length() <= N_STRAND) {
                    ss = ss.replace('B', '-');
                }

                builder.append(ss);
                i = l - 1;
            } else if (c == 'H') {

                StringBuilder sBuilder = new StringBuilder();
                int l = i;
                char cH = 'H';
                do {

                    sBuilder.append(cH);
                    l++;
                    if (l < nres)
                        cH = secondaryStructure.charAt(l);

                } while (cH == 'H' && l < nres);

                String ss = sBuilder.toString();

                if (ss.length() <= N_HELIX) {
                    ss = ss.replace('H', '-');
                }

                builder.append(ss);
                i = l - 1;

            }

        }

        return builder.toString();
    }

    /**
     * Secondary structure generated from PDB be notice that we DO NOT use
     * information from DSSP file incase of unparseable DSSP file
     *
     * @return secondary structure in sequence
     */
    public String getCombinedSS(Atom[] caAtoms, String chain) {
//        StringBuilder builder = new StringBuilder();
        String strSS = "";
        String strSSPDB = this.getSS(caAtoms);
        try {
            strSS = this.getSS(chain, caAtoms);
        } catch (Exception ex) {
            // do something here
            ex.printStackTrace();
        }

        if (strSS.length() > 0)
            return strSS;//
            // in the case we cannot parse the DSSP file, we use the information from PDB
        else return strSSPDB;


        // StrucState[] strucStates = this.getStrucState(chain, caAtoms);
        // Sequence3D sequence3D = new Sequence3D();
        // String seq3d = "";
        // try {
        // String strAccs = "";
        // for (int i = 0; i < caAtoms.length; i++) {
        //
        // strAccs += "1";
        // }
        //
        // seq3d = sequence3D.getSA(strucStates, strAccs);
        //
        // } catch (Exception e) {
        // // TODO Auto-generated catch block
        // e.printStackTrace();
        // }
        // seq3d = seq3d.toUpperCase();
        // System.out.println(seq3d);
        // combined all.

//        for (int i = 0; i < caAtoms.length; i++) {
//
//            if (strSSPDB.charAt(i) == 'B' || strSSDSSP.charAt(i) == 'B')
//                builder.append('B');
//            else if (strSSPDB.charAt(i) == 'H' || strSSDSSP.charAt(i) == 'H')
//                builder.append('H');
//            else
//                builder.append("-");
//        }

//        return builder.toString();
    }

    /**
     * Secondary structure generated from PDB be notice that we DO NOT use
     * information from DSSP file
     *
     * @return secondary structure in sequence
     */
    public String getSS(Atom[] caAtoms) {

        StringBuilder builder = new StringBuilder();

        for (int i = 0; i < caAtoms.length; i++) {
            try {
                SecStrucInfo sec = (SecStrucInfo) (caAtoms[i].getGroup())
                        .getProperty(Group.SEC_STRUC);

                String ss = sec.getType().toString().trim();


                if (ss.equals("H") || ss.equals("G") || ss.equals("I"))
                    builder.append('H');
                else if (ss.equals("E") || ss.equals("B"))
                    builder.append('B');
                else
                    builder.append('-');

//                if(sec.get)
//                } else {
//                    builder.append("-");
//                }

            } catch (Exception ex) {
                // NOTE : if this protein does not have any information about
                // secondary structure

                builder.append("-");
            }

        }


//        for (int i = 0; i < caAtoms.length; i++) {
//            try {
//                Map<String, String> sec = ((AminoAcid) caAtoms[i].getGroup())
//                        .getSecStruc();
//                if (sec.containsKey("PDB_AUTHOR_ASSIGNMENT")) {
//                    String str2nd = sec.get("PDB_AUTHOR_ASSIGNMENT");
//                    if (str2nd.equals("STRAND")) {
//                        builder.append("B");
//                    } else if (str2nd.equals("HELIX")) {
//                        builder.append("H");
//                    }
//
//                } else {
//                    builder.append("-");
//                }
//
//            } catch (Exception ex) {
//                // NOTE : if this protein does not have any information about
//                // secondary structure
//
//                builder.append("-");
//            }
//
//        }
        return builder.toString();
    }

    /**
     * @param chainName
     * @return
     */
    public String getSS(String chainName, Atom[] atoms) {

        StringBuilder builder = new StringBuilder();
        Map<Integer, Character> mapAtomDSSP = new HashMap<Integer, Character>();

        for (DSSPDataLine line : c) {
            if (chainName.equals(line.chainLetter + "")) {
                char ss = line.structure;
                mapAtomDSSP.put(line.seqPosition, ss);
            }
        }

        for (Atom atom : atoms) {
            char ss = '-';
            if (mapAtomDSSP.containsKey(atom.getGroup().getResidueNumber()
                    .getSeqNum()))
                ss = mapAtomDSSP.get(atom.getGroup().getResidueNumber()
                        .getSeqNum());
            if (ss == 'H' || ss == 'G' || ss == 'I')
                builder.append('H');
            else if (ss == 'E' || ss == 'B')
                builder.append('B');
            else
                builder.append('-');

        }

        return builder.toString();
    }

    public String getRelativeACC(String chainName, Atom[] atoms) {

        Map<String, Double> map = RelativeASA.getAASurface();

        StringBuilder builder = new StringBuilder();

        Map<Integer, Integer> mapAtomDSSP = new HashMap<Integer, Integer>();

        for (DSSPDataLine line : c) {
            if (chainName.equals(line.chainLetter + "")) {

                double acc = 1.0;
                if (map.containsKey("" + line.residueLetter))
                    acc = (double) line.acc / map.get("" + line.residueLetter);
                // System.out.println(index + ":" + acc);
                // index++;
                if (acc <= THRES_ACC)
                    mapAtomDSSP.put(line.seqPosition, 0);
                else
                    mapAtomDSSP.put(line.seqPosition, 1);
            }
        }

        for (Atom atom : atoms) {
            int ss = 1;
            if (mapAtomDSSP.containsKey(atom.getGroup().getResidueNumber()
                    .getSeqNum()))
                ss = mapAtomDSSP.get(atom.getGroup().getResidueNumber()
                        .getSeqNum());

            builder.append(ss);

        }

        return builder.toString();
    }

    public StrucState[] getStrucState(String chainName, Atom[] atoms) {

        StrucState[] pr = new StrucState[atoms.length];
        for (int i = 0; i < atoms.length; i++) {
            pr[i] = new StrucState();
        }

        Map<Integer, StrucState> mapAtomDSSP = new HashMap<Integer, StrucState>();

        for (DSSPDataLine line : c) {
            if (chainName.equals(line.chainLetter + "")) {
                StrucState struState = new StrucState();
                struState.setPhi(line.phi);
                struState.setPsi(line.psi);
                struState.setAlpha(line.alpha);
                struState.setKappa(line.kappa);

                mapAtomDSSP.put(line.seqPosition, struState);

            }
        }

        for (int i = 0; i < atoms.length; i++) {
            if (mapAtomDSSP.containsKey(atoms[i].getGroup().getResidueNumber()
                    .getSeqNum()))
                pr[i] = mapAtomDSSP.get(atoms[i].getGroup().getResidueNumber()
                        .getSeqNum());

        }
        return pr;
    }

    // public Integer[] getACCValues(String chainName) {
    // List<Integer> list = new ArrayList<Integer>();
    // for (DSSPDataLine line : c) {
    // if (chainName.equals(line.chainLetter + ""))
    // list.add(line.acc);
    // }
    // return (Integer[]) (list.toArray(new Integer[list.size()]));
    // }

    private void parser(String paramString, boolean paramBoolean) {
        for (; ; ) {
            if (this.d < 2) {
                // if (READ_CLASSIC_FORMAT) {
                // e.fine("Attempting CLASSIC file format");
                // } else {
                // e.fine("Attempting NEW file format");
                // }
            } else {
                break;
                // e.severe("No more file format variants to try - is it really a DSSP file?");
                // System.exit(0);
            }
            this.d = ((byte) (this.d + 1));
            String str = "";
            try {
                BufferedReader localBufferedReader = new BufferedReader(
                        new FileReader(paramString));
                do {
                    if ((str = localBufferedReader.readLine()) == null) {
                        // e.severe("The file seems to be empty or the format is incorrect: "
                        // + paramString);
                        // System.exit(0);
                    }
                } while (str.charAt(2) != '#');
                while ((str = localBufferedReader.readLine()) != null) {
                    // System.out.println(str);
                    if ((str.length() >= 120) && (str.charAt(13) != '!')) {
                        DSSPDataLine localDSSPDataLine = new DSSPDataLine(str,
                                paramBoolean);

                        this.c.add(localDSSPDataLine);
                        if (localDSSPDataLine.position > this.a) {
                            this.a = localDSSPDataLine.position;
                        }
                    }
                }
                localBufferedReader.close();
                return;
            } catch (IOException localIOException) {
                logger.error("error", localIOException);


                // System.exit(0);
            } catch (Exception localException) {
                //  localException.printStackTrace();
                logger.error(localException.toString());
            }
            paramBoolean = !paramBoolean;
            paramString = paramString;

        }

    }

}