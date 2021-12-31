package org.cnrs.crbm.lib.utils;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.repeats.RepeatFinder;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 10/16/2015.
 */
public class FastaUtils {


    /**
     * File format
     * 4qrl_A
     * 2bol_A
     *
     * @param fileIn
     * @param fileOut
     */
    public static void exportFastaFileFromListPDB(String fileIn, String fileOut) throws IOException {

        // read input
//        PrintWriter writer = new PrintWriter("C:/Users/pdoviet/Desktop/TAPO_NoTRs.fasta");
        PrintWriter writer = new PrintWriter(fileOut);
        //List<String> pdbs = DataIO.readLines("data/resultLargeScale/TaPoTRs.txt");
        List<String> pdbs = DataIO.readLines(fileIn);
        Row row = null;
        LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
                .readFastaProteinSequence(new File("data/nrdb/pdb_seqres.txt"));
        for (Map.Entry<String, ProteinSequence> entry : a.entrySet()) {
            // System.out.println(entry.getValue().getOriginalHeader() + "="
            // + entry.getValue().getSequenceAsString());
//            if(pdbs.contains(entry.getKey())){
//                System.out.println(entry.getKey());
//            }else
//                System.err.println(entry.getKey());

            row = new Row();
            row.setLength(entry.getValue().getLength());
            row.setProtein(entry.getValue().getOriginalHeader().substring(0, 6));
            String pdb = row.getPdbCode().toLowerCase() + "_" + row.getPdbChain();
            if (pdbs.contains(pdb)) {
                //System.out.println(pdb);
                //System.out.println(entry.getValue().getSequenceAsString());
                writer.write(">" + pdb + "\n");
                writer.write(entry.getValue().getSequenceAsString() + "\n");
            }

        }// end

        writer.close();

    }

    /**
     * File format
     * pdbId	regions	predLen	predRU	classId	cathId	ranges	topologyID
     * 4qrl_A.4_107	30-133	25.5	4	2	unk	20-40	unk
     * 2bol_A.256_300	258-314	21	2	2	unk	20-40	unk
     * 3am2_A.194_270	230-306	26	3	2	unk	20-40	unk
     *
     * @param fileIn
     * @param fileOut
     */
    public static void exportFastaFileFromListPDBRegion(String fileIn, String fileOut) {
        List<String> pdbs = DataIO.readLines(fileIn);
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = pdbs.size();
        StringBuffer buffer = new StringBuffer();
        for (String pdb : pdbs) {
            bar.update(process, sizeOfProcess);
            process++;
            pdb = pdb.split("\t")[0];
            String pdbID = pdb.split("\\.")[0];
            String regions = pdb.split("\\.")[1];
            int start = Integer.parseInt(regions.split("_")[0]);
            int end = Integer.parseInt(regions.split("_")[1]);
            String pdbCode = pdbID.substring(0, 4);
            String pdbChain = pdbID.substring(5, 6);
            RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
            buffer.append(">" + pdb + "\n");
            buffer.append(finder.getStrSeq().substring(start, end + 1) + "\n");

        }

        DataIO.writeToFile(buffer.toString(), fileOut);
    }
}
