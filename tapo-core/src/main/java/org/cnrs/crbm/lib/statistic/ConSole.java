package org.cnrs.crbm.lib.statistic;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.parser.Parser;
import org.jsoup.select.Elements;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by pdoviet on 2/5/2015.
 */
public class ConSole {


    private static String WEBSERVER = "http://console.sanfordburnham.org.";

    public static void main(String[] args) throws Exception {

        ConSole conSole = new ConSole();
        conSole.runEvaluation();

        //Row row = new Row();
        //row.setProtein("1ifp_A");
        //conSole.submitQuery(row);
        //conSole.convertXmltoFile(row);
        //conSole.downloadResultFromServer(row, "res1lxaA20150205003241.xml");


    }

    private void convertXmltoFile(Row row) throws IOException {

        String xml = DataIO.readFile("output/console/" + row.getProtein() + ".xml");
        Document xmlDoc = Jsoup.parse(xml, "", Parser.xmlParser());
        Elements pDBClassification = xmlDoc.select("PDBClassification");

        //System.out.println(pDBClassification.attr("MeanDistance"));

        String meanDistance = pDBClassification.attr("MeanDistance");
        String numberSampledResidues = xmlDoc.select("Summary").attr("NumberSampledResidues");
        String numberSolenoids = xmlDoc.select("Summary").attr("NumberSolenoids");

        String structureClassification = "[" + xmlDoc.select("StructureClassification").attr("Criterion0.33");
        structureClassification += " " + xmlDoc.select("StructureClassification").attr("Criterion0.5");
        structureClassification += " " + xmlDoc.select("StructureClassification").attr("Criterion0.66");
        structureClassification += " " + xmlDoc.select("StructureClassification").attr("isRing") + "]";


        String units = "";
        String unitsFasta = xmlDoc.select("UnitsFasta").text();
        //System.out.println(unitsFasta);
        for (Element e : xmlDoc.select("Unit")) {
            units += e.attr("StartResID") + "-" + e.attr("EndResID") + ";";
        }

        String content = row.getProtein() + "\t" + meanDistance + "\t" + numberSampledResidues + "\t" + numberSolenoids + "\t" + structureClassification + "\t" + units + "\t" + unitsFasta;
        DataIO.writeToFile(content, "output/console/converted/" + row.getProtein() + ".tab");

    }


    public void runEvaluation() throws Exception {

        // read pdbList

        File folder = new File("output/console/");

        File[] listOfFiles = folder.listFiles();
        Set<String> set = new HashSet<String>();

        for (File file : listOfFiles) {
            if (file.isFile()) {
                set.add(file.getName());
            }
        }

        set.add("3ikm_A.xml");
        set.add("3dtp_B.xml");
        set.add("1ea0_A.xml");
        set.add("1ofe_B.xml");
        set.add("2b39_A.xml");


        List<Row> rows = DataIO.getListProteinFromFile("data/pdbList.in");
        for (Row row : rows) {
            //System.out.println(row.getProtein());
            if (set.contains(row.getProtein() + ".xml"))
                continue;

            RepeatFinder finder = new RepeatFinder(row.getPdbCode(), row.getPdbChain());
//            if (finder.getAtoms().length > 200)
//                continue;
            this.submitQuery(row);
            this.convertXmltoFile(row);
            System.out.printf("\r" + row.getProtein() + " is completed!");
        }

    }

    public void submitQuery(Row row) throws Exception {

        System.out.printf("\r" + row.getProtein() + " is submitting...");
        String queryUrl = WEBSERVER + "/uploadPDBFile.php?pdbID=" + row.getPdbCode() + "&chain=" + row.getPdbChain();

        URL url = new URL(queryUrl);
        URLConnection con = url.openConnection();

        BufferedReader rd = new BufferedReader(new InputStreamReader(
                con.getInputStream()));
        StringBuffer sb = new StringBuffer();
        String line;
        while ((line = rd.readLine()) != null) {
            sb.append(line + "\n");
        }

        rd.close();
        Document doc = Jsoup.parse(sb.toString());
        //Document doc = Jsoup.connect(queryUrl).get();
        Element link = doc.select("a").first();
        String linkHref = link.attr("href");
        this.downloadResultFromServer(row, linkHref);

    }


    public void downloadResultFromServer(Row row, String fileName) {

        String fileDir = "output/console/" + row.getProtein() + ".xml";
        // check file exist
        if (!new File(fileDir).exists()) {
            // load
            try {
                String pdbfasta = WEBSERVER + "/" + fileName;
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

            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }


    }

}
