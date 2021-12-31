package apps;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureImpl;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class RepeatDBAnalysis {

    public static void main(String[] args) throws Exception {
//        new RepeatDBAnalysis().saveDataToMaxCluster();

        // new RepeatDBAnalysis().proteinClassification();
        // new RepeatDBAnalysis().run();
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");

        List<Row> rows = new ArrayList<Row>();
        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            //if (r.getAnnlevel().equals("Detailed")) {
            System.out.println(r.getPdbCode() + "_" + r.getPdbChain());
            //}
        }

//        new RepeatDBAnalysis().compareTAPOvsRepeatsDB();

    }


    public void saveDataToMaxCluster() throws Exception {

        List<String> lines = DataIO.readLines("data/pdbList.in");

        List<Row> rows = new ArrayList<Row>();
        for (String line : lines) {

            if (line.startsWith("#"))
                continue;

            //System.out.println(line.substring(0,5));
            Row row = new Row();
            row.setProtein(line.substring(0, 6));
            rows.add(row);

        }
        this.saveToTmp(rows);


    }


    public void compareTAPOvsRepeatsDB() throws Exception {

        Set<String> set = new HashSet<String>();
        List<Row> rows = this
                .getListProteinFromLargeScaleAnalysis("data/tapo/tapo_pdb_July_1_2011.o");
        for (Row row : rows) {
            //System.out.println(row.getProtein());
            set.add(row.getProtein());
        }


        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");


        Set<String> set2 = new HashSet<String>();
        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            String proteinCode = r.getPdbCode() + "_" + r.getPdbChain();

            if (!set.contains(proteinCode) && !set2.contains(proteinCode)) {
                System.out.println(proteinCode + "\t" + r.getAnnlevel() + "\t" + r.getStrclass());
            }
            set2.add(proteinCode);

        }


    }

//    public void proteinClassification() throws Exception {
//
//        ProClassify proClassify = new ProClassify();
//        List<String> lines = DataIO.readLines("data/pdbClassification.o");
//        Map<String, String> map = this.getListProteinOutput("data/finalPDB.o");
//        for (String line : lines) {
//            try {
//                String pdb = line.substring(0, 6);
//                String pdbCode = pdb.substring(0, 4);
//                String pdbChain = pdb.substring(5, 6);
//                String content = map.get(pdb);
//
//                System.out.println(line
//                        + "\t"
//                        + content
//                        + "\t"
//                        + proClassify.getRepresentativeDomains(pdbCode,
//                        pdbChain));
//
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }
//
//    }

    public void run() throws Exception {

        List<Row> repeatDBs = this.getRepeatDB();
        ProteinCSVReader csvReader = new ProteinCSVReader();

        List<Row> rowsFrancois = csvReader.getData("data/RepeatDatalastest.csv");

        Set<String> pdbsTRs = new HashSet<String>();
        List<Row> rows = this
                .getListProteinFromLargeScaleAnalysis("data/finalPDB.o");

        // List<Row> nrpdb2013 = DataIO
        // .getListProteinFromFile("data/nrpdb2013.txt");
        // System.out.println("nrpdb2013 : " + nrpdb2013.size());
        System.out.println("RepeatDB : " + repeatDBs.size());
        System.out.println("Seqs Method TRs : " + rowsFrancois.size());
        System.out.println("3D Method TRs : " + rows.size());

        for (Row row : rows) {
            pdbsTRs.add(row.getProtein());
        }

        System.out.println();
        // check overlapp with Francois
        int count = 0;
        for (Row row : rowsFrancois) {
            //System.out.println(row.getProtein());

            if (pdbsTRs.contains(row.getProtein())) {
                count++;
                pdbsTRs.remove(row.getProtein());

            }
        }

        System.exit(0);
        System.out.println("Francois overlapping: " + count + "/"
                + rowsFrancois.size());

        count = 0;
        for (Row row : repeatDBs) {
            System.out.println(row);
            if (pdbsTRs.contains(row.getProtein())) {
                count++;
                pdbsTRs.remove(row.getProtein());

            }
        }

        System.out.println("RepeatDB overlapping: " + count + "/"
                + repeatDBs.size());

        System.out.println("new pdb TRs: " + pdbsTRs.size());

        StringBuffer buffer = new StringBuffer();
        for (String pdb : pdbsTRs) {
            // System.out.println(pdb);
            buffer.append(pdb + "\n");

        }

        // DataIO.writeToFile(buffer.toString(), "data/pdbList.in");

        // System.out.println("Resutls: " + a);
        // saveToTmp(saveRows);

    }

    public void saveToTmp(List<Row> rows) throws StructureException {

        for (Row row : rows) {
            // if (!pdbsTRs.contains(row.getProtein()))
            // System.out.println(row.getProtein());
            //

            Structure structure = PdbTools.getStructureFromLocalPdb(row
                    .getPdbCode());

            Chain c = structure.getChainByPDB(row.getPdbChain());

            Structure newStruct = new StructureImpl();
            newStruct.addChain(c);

            DataIO.writeToFile(newStruct.toPDB(),
                    "F:\\Cygwin\\home\\pdoviet\\BioApps\\MaxCluster\\pdbs\\"
                            + row.getProtein() + ".pdb");

        }

    }

    public Map<String, String> getListProteinOutput(String fileDir) {

        Map<String, String> map = new HashMap<String, String>();
        List<String> lines = DataIO.readLines(fileDir);
        for (String line : lines) {
            int nrTRs = Integer.parseInt(line.split("\t")[1]);
            String annotation = line.split("\t")[2];
            int len = Integer.parseInt(line.split("\t")[3]);
            if (nrTRs >= 2)
                map.put(line.substring(0, 6), nrTRs + "\t" + annotation + "\t"
                        + len);
        }

        return map;

    }

    public List<Row> getListProteinFromLargeScaleAnalysis(String fileDir) {

        List<Row> rows = new ArrayList<Row>();

        List<String> lines = DataIO.readLines(fileDir);
        for (String line : lines) {
            if (line.length() > 0) {
                Row row = new Row();
                row.setProtein(line.substring(0, 6));
//			int nrTRs = Integer.parseInt(line.split("\t")[1]);
//			if (nrTRs == 2)
                rows.add(row);
            }
        }

        return rows;

    }

    public List<Row> getRepeatDB() {

        List<Row> rows = new ArrayList<Row>();

        BufferedReader br = null;

        try {

            String sCurrentLine;

            br = new BufferedReader(new FileReader("data/RDB-dataset.tab"));

            br.readLine();
            while ((sCurrentLine = br.readLine()) != null) {
                // System.out.println(sCurrentLine);
                Row row = new Row();
                row.setProtein(sCurrentLine.substring(0, 4) + "_"
                        + sCurrentLine.substring(4, 5));

                rows.add(row);

            }

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null)
                    br.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

        return rows;

    }

}
