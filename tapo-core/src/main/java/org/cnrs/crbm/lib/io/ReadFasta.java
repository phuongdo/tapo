package org.cnrs.crbm.lib.io;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;
import java.util.Map.Entry;


public class ReadFasta {


    public static void main(String[] args) throws Exception {
//        ReadFasta fasta = new ReadFasta();
//        Map<String, String> map = fasta.getFastaFileInSpecialCase("data/pdbs.in.ali");
//
//
//        System.out.println(map.entrySet().iterator().next().getValue());
//        for (Map.Entry<String, String> entry : map.entrySet()) {
//            System.out.println(entry.getKey() + "/" + entry.getValue());
//        }

        ReadFasta readFasta = new ReadFasta();
        readFasta.readFastaTapoFormatTest();
    }


    /**
     * TESTING METHOD
     *
     * @throws Exception
     */
    public void readFastaTapoFormatTest() throws Exception {


        Map<String, TaPoFastaFormat> fastaTRs = this.readTaPoFastaFormat("data/rossmann/rossmann.TRs.o");
        TaPoFastaFormat value = fastaTRs.get("2e0p_A");
        System.out.println(value);

    }

    /**
     * @return the alphabet string of one pdb structure
     * @throws Exception
     */
    public Set<String> getAllPdbs(String dir) throws Exception {
        // Try with the FastaReaderHelper
        Set<String> pdbs = new HashSet<String>();
        LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
                .readFastaProteinSequence(new File(dir));
        for (Entry<String, ProteinSequence> entry : a.entrySet()) {
            // System.out.println(entry.getValue().getOriginalHeader() + "="
            // + entry.getValue().getSequenceAsString());
            pdbs.add(entry.getValue().getOriginalHeader().substring(0, 4));

        }
        return pdbs;

    }

    /**
     * @return the alphabet string of one pdb structure
     * @throws Exception
     */
    public Set<String> getAllPdbsWithChain(String dir) throws Exception {
        // Try with the FastaReaderHelper
        Set<String> pdbs = new HashSet<String>();
        LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
                .readFastaProteinSequence(new File(dir));
        for (Entry<String, ProteinSequence> entry : a.entrySet()) {
            // System.out.println(entry.getValue().getOriginalHeader() + "="
            // + entry.getValue().getSequenceAsString());
            pdbs.add(entry.getValue().getOriginalHeader().substring(0, 6));

        }
        return pdbs;

    }

    /**
     * @return the alphabet string of one pdb structure
     * @throws Exception
     */
    public List<Row> getRowAllPdbsWithChain(String dir) throws Exception {
        // Try with the FastaReaderHelper
        List<Row> pdbs = new ArrayList<Row>();
        // int i = 0;
        LinkedHashMap<String, ProteinSequence> a = FastaReaderHelper
                .readFastaProteinSequence(new File(dir));
        for (Entry<String, ProteinSequence> entry : a.entrySet()) {
            Row row = new Row();
            row.setLength(entry.getValue().getLength());
            row.setProtein(entry.getValue().getOriginalHeader().substring(0, 6));
            pdbs.add(row);
            // i++;
            // if (i > 30)
            // break;

        }
        return pdbs;

    }


    public Map<String, String> getFastaFileInSpecialCase(String fileDir) throws Exception {
        // Try with the FastaReaderHelper

        Map<String, String> ss = new HashMap<String, String>();
        BufferedReader br = new BufferedReader(new FileReader(fileDir));

        List<String> name = new ArrayList<String>();
        List<String> details = new ArrayList<String>();
        String line;
        // Now read lines of text: the BufferedReader puts them in lines,
        // the InputStreamReader does Unicode conversion, and the
        // GZipInputStream "gunzip"s the data from the FileInputStream.
        StringBuilder builder = new StringBuilder();
        while ((line = br.readLine()) != null) {


            if (line.startsWith(">")) {
                // System.out.println(line);
                name.add(line);

            } else {
                builder = new StringBuilder();
                builder.append(line);
                while ((line = br.readLine()) != null) {

                    if (line.startsWith(">")) {

                        details.add(builder.toString().replaceAll("\0", "-"));
                        name.add(line);
                        //name.add(line.replace(">", ""));
                        break;
                    } else
                        builder.append(line);
                }

            }

        }
        details.add(builder.toString());

        for (int i = 0; i < name.size(); i++) {

            String aname = name.get(i).replace(">", "");
            String adescription = details.get(i);


            ss.put(aname, adescription.toString());
        }

        br.close();
        return ss;
    }

    /**
     * Read tapo output based on FASTA format.
     *
     * @param fileDir
     * @return
     * @throws Exception
     */
    public Map<String, TaPoFastaFormat> readTaPoFastaFormat(String fileDir) throws Exception {
        // Try with the FastaReaderHelper

        Map<String, TaPoFastaFormat> ss = new HashMap<String, TaPoFastaFormat>();
        BufferedReader br = new BufferedReader(new FileReader(fileDir));

        List<String> name = new ArrayList<String>();
        List<String> details = new ArrayList<String>();
        String line;
        // Now read lines of text: the BufferedReader puts them in lines,
        // the InputStreamReader does Unicode conversion, and the
        // GZipInputStream "gunzip"s the data from the FileInputStream.
        StringBuilder builder = new StringBuilder();
        while ((line = br.readLine()) != null) {


            if (line.startsWith(">")) {
                // System.out.println(line);
                name.add(line);

            } else {
                builder = new StringBuilder();
                builder.append(line + "|");
                while ((line = br.readLine()) != null) {
                    if (line.startsWith(">")) {
                        details.add(builder.toString().replaceAll("\0", "-"));
                        name.add(line);
                        //name.add(line.replace(">", ""));
                        break;
                    } else
                        builder.append(line + "|");
                }

            }

        }
        details.add(builder.toString());

        for (int i = 0; i < name.size(); i++) {

            String aheader = name.get(i).replace(">", "");
            String adescription = details.get(i);
            TaPoFastaFormat taPoFastaFormat = new TaPoFastaFormat(aheader, adescription);
            ss.put(taPoFastaFormat.getPdbID(), taPoFastaFormat);
            //taPoFastaFormat.setHeader(aheader);
//            ss.put(aname, adescription.toString());
        }

        br.close();
        return ss;
    }

}
