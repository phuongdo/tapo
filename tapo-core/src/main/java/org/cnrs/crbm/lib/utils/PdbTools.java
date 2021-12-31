package org.cnrs.crbm.lib.utils;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;


import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.LocalPDBReader;

public class PdbTools {

    private final static String PDB_SERVICE = "http://www.pdb.org/pdb/files";
    private final static String FASTA_SERVICE = "http://www.rcsb.org/pdb/files/fasta.txt?structureIdList";// $pdbid

    @Deprecated
    public static String downloadPDB(String pdbCode) {

        // check file exist

        String fileDir = Dir.DB_TMP_DIR + "/" + pdbCode + ".pdb";
        if (!new File(fileDir).exists()) {
            // load
            try {
                String pdbGZ = PDB_SERVICE + "/" + pdbCode + ".pdb.gz";
                URL url = new URL(pdbGZ);
                URLConnection con = url.openConnection();
                GZIPInputStream gis = new GZIPInputStream(con.getInputStream());
                FileOutputStream fos = new FileOutputStream(fileDir);
                byte[] buffer = new byte[1024];
                int len;
                while ((len = gis.read(buffer)) != -1) {
                    fos.write(buffer, 0, len);
                }
                // close resources

                fos.close();
                gis.close();

            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        // remove this
        // copyFileToJmolWebApplication(pdbCode);

        return fileDir;

    }

    public static String copyFileToJmolWebApplication(String pdbCode) {
        // check file exist

        String fileDir = Dir.JMOL_WEB_DATA + "/" + pdbCode + ".pdb.gz";
        if (!new File(fileDir).exists()) {
            // load
            try {
                String pdbGZ = PDB_SERVICE + "/" + pdbCode + ".pdb.gz";
                URL url = new URL(pdbGZ);
                URLConnection con = url.openConnection();

                BufferedInputStream bis = new BufferedInputStream(
                        con.getInputStream());
                FileOutputStream fos = new FileOutputStream(fileDir);

                // GZIPInputStream gis = new
                // GZIPInputStream(con.getInputStream());

                byte[] buffer = new byte[1024];
                int len;
                while ((len = bis.read(buffer)) != -1) {
                    fos.write(buffer, 0, len);
                }
                // close resources

                fos.close();
                bis.close();

            } catch (Exception e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        return fileDir;

    }

    public static String insertToFastaDB(Structure structure) {
        String pdbCode = structure.getPDBCode().toLowerCase();
        String fileDir = Dir.FASTA_LOCAL + "/" + pdbCode.substring(1, 3) + "/"
                + pdbCode + ".fasta";
        // check file exist
        if (!new File(fileDir).exists()) {
            // load
            try {

                String dir = Dir.FASTA_LOCAL + "/" + pdbCode.substring(1, 3);
                File theDir = new File(dir);
                if (!theDir.exists()) {
                    // make a directory before save file
                    theDir.mkdir();
                }

                StringBuilder sb = new StringBuilder();

                for (Chain c : structure.getChains()) {
                    sb.append(">" + pdbCode.toUpperCase() + ":"
                            + c.getChainID() + "|PDBID|CHAIN|ATOM SEQUENCE\n");
                    // print biological sequence
                    sb.append(c.getAtomSequence() + "\n");
                }

                DataIO.writeToFile(sb.toString(), fileDir);

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        return fileDir;

    }

    public static String downloadFastaFromServer(String pdbCode) {

        String fileDir = Dir.DB_TMP_DIR + "/" + pdbCode + ".fasta";
        // check file exist
        if (!new File(fileDir).exists()) {
            // load
            try {
                String pdbfasta = FASTA_SERVICE + "=" + pdbCode;
                URL url = new URL(pdbfasta);
                URLConnection con = url.openConnection();

                BufferedReader rd = new BufferedReader(new InputStreamReader(
                        con.getInputStream()));
                StringBuffer sb = new StringBuffer();
                String line;
                while ((line = rd.readLine()) != null) {
                    sb.append(line + "\n");
                }

                rd.close();
                DataIO.writeToFile(sb.toString(), fileDir);

                // DataInputStream gis = new
                // DataInputStream(con.getInputStream());
                // FileOutputStream fos = new FileOutputStream(fileDir);
                // byte[] buffer = new byte[1024];
                // int len;
                // while ((len = gis.read(buffer)) != -1) {
                // fos.write(buffer, 0, len);
                // }
                // // close resources
                //
                // fos.close();
                // gis.close();

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        return fileDir;

    }

    public static Structure getStructureFromFile(String pdbFile) {
        FileParsingParameters params = new FileParsingParameters();
        params.setParseSecStruc(true);
        PDBFileReader pdbreader = new PDBFileReader();
        pdbreader.setFileParsingParameters(params);
        Structure struc1 = null;
        try {
            struc1 = pdbreader.getStructure(pdbFile);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return struc1;
    }

    /**
     * a missing PDB id be fetched automatically from the FTP servers
     *
     * @param pdbCode
     * @return PDB structure
     */
    public static Structure getStructureFromLocalPdb(String pdbCode) {
        LocalPDBReader localPDBReader = new LocalPDBReader();
        return localPDBReader.getPdb(pdbCode);
    }


    /**
     * Load a structure from file
     *
     * @param filedDir
     * @return
     */
    public static Structure loadStructureFromFile(String filedDir) {
        LocalPDBReader localPDBReader = new LocalPDBReader();
        return localPDBReader.loadStructureFromFile(filedDir);
    }

    /**
     * only amino acid
     *
     * @param
     * @return
     */
    public static Atom[] getAtomCAArray(Chain chain) {
        Atom[] caAtoms = StructureTools.getAtomCAArray(chain);
        List<Atom> listAtoms = new ArrayList<Atom>();
        for (Atom atom : caAtoms) {
            if (atom.getGroup() instanceof AminoAcid) {
                listAtoms.add(atom);
            }
        }
        return listAtoms.toArray(new Atom[listAtoms.size()]);
    }

    public static Atom[] getAtomCBArray(Atom[] caAtoms) {
        int i = 0;
        Atom[] arrAtomCp = new Atom[caAtoms.length];
        for (Atom atom : caAtoms) {
            arrAtomCp[i] = atom;
            if (atom.getGroup() instanceof AminoAcid) {
                AminoAcid aa = (AminoAcid) atom.getGroup();
                for (Atom a : aa.getAtoms()) {
                    if (a.getName().equals("CB")) {
                        arrAtomCp[i] = a;
                    }
                }
            }

            i++;
        }
        return arrAtomCp;
    }

    /**
     * Write a pdb file from structure object.
     */

    public static void writePDB(Atom[] atoms, String pdbCode, String pdbChain,
                                int start, int end, String outputfile) {
        try {

            Structure struc = new StructureImpl();
            Chain c = new ChainImpl();

            for (int i = start; i <= end; i++) {

                try {
                    Group group = atoms[i].getGroup();
                    c.addGroup(group);
                } catch (Exception ex) {
                    // Not show exception here
                }
            }

            struc.addChain(c);

            // String dir = "C:/Users/pdoviet/Desktop";
            // // write to pdb file
            // String outputfile = dir + "/" + pdbCode.toLowerCase() + pdbChain
            // + "_" + start + "_" + end + ".pdb";
            FileOutputStream out = new FileOutputStream(outputfile);
            PrintStream p = new PrintStream(out);
            p.println(struc.toPDB());
            p.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    public static void writePDB(List<Atom> atoms, String outputfile) {
        try {

            Structure struc = new StructureImpl();
            Chain c = new ChainImpl();

            StringBuffer buffer = new StringBuffer();
            for (int i = 0; i < atoms.size(); i++) {

                try {
                    Group group = atoms.get(i).getGroup();
                    atoms.get(i).setPDBserial(i + 1);
                    buffer.append(atoms.get(i).toPDB() + "\n");
                    c.addGroup(group);
                } catch (Exception ex) {
                    // Not show exception here
                }
            }

            struc.addChain(c);
            FileOutputStream out = new FileOutputStream(outputfile);
            PrintStream p = new PrintStream(out);
            //p.println(struc.toPDB());
            p.println(buffer.toString());
            p.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    public static String getSeqResSequence(Atom[] atoms) {
        StringBuffer buffer = new StringBuffer();
        for (Atom atom : atoms) {
            buffer.append(atom.getGroup().getChemComp().getOne_letter_code().charAt(0));
        }
        return buffer.toString();
    }

    public static int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {
            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;
        }
        return 0;
    }

    public static int getResSeq(Atom[] atoms, int pos) {
        return atoms[pos].getGroup().getResidueNumber().getSeqNum();
    }

}
