package org.cnrs.crbm.lib.statistic;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.networking.MultipartUtility;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * Created by pdoviet on 5/29/2015.
 */
public class ProStripServer {

    public static void main(String[] args) throws Exception {

        ProStripServer proStripServer = new ProStripServer();

        try {
            //proStripServer.runEvaluation();
            proStripServer.exportResult();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void exportResult() throws Exception{


        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");

//        File folder = new File("output/prostrip/");
//        File[] listOfFiles = folder.listFiles();
//        Set<String> set = new HashSet<String>();
//        for (File file : listOfFiles) {
//            if (file.isFile()) {
//                set.add(file.getName());
//
//
//            }
//        }

        for (String row : rows) {
            if (row.startsWith("#"))
                continue;
            String[] datas = row.split("\t");

            if(datas[1].equals("0")) {

                System.out.print(DataIO.readFile("output/prostrip/"+datas[0]+".pro"));

            }

        }





    }

    public void runEvaluation() throws Exception {

        List<Row> rows = DataIO.getListProteinFromFile("data/pdbList.in");
        Random random = new Random();
        File folder = new File("output/prostrip/");
        File[] listOfFiles = folder.listFiles();
        Set<String> set = new HashSet<String>();
        for (File file : listOfFiles) {
            if (file.isFile()) {
                set.add(file.getName());
            }
        }

        for (Row row : rows) {
            try {
                String result = this.submitQuery(row);

                String filePathString = "output/prostrip/" + row.getProtein() + ".pro";
                File f = new File(filePathString);
                if (set.contains(row.getProtein() + ".pro"))
                    continue;

                String urlResults = this.submitQuery(row);

                System.out.print(".");

                if (urlResults.equals(""))
                    break;

                if (f.exists() && !f.isDirectory()) {
                    continue;
                } else
                    DataIO.writeToFile(result, "output/prostrip/" + row.getProtein() + ".pro");

//                long minutes = 100 * 1000;
//                long randomSecond = random.nextInt(30) * 1000;
//                Thread.sleep(minutes + randomSecond);
            } catch (Exception ex) {

            }
        }

    }

    private String submitQuery(Row row) throws Exception {

        String charset = "UTF-8";
        String requestURL = "http://cluster.physics.iisc.ernet.in/cgi-bin/prostrip/frame.pl";

        //StringBuffer stringBuffer = new StringBuffer();
        String results = "";
        String urlResult = "";
        try {
            MultipartUtility multipart = new MultipartUtility(requestURL, charset, "ProStrip");

            multipart.addHeaderField("User-Agent", "Agent:Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/43.0.2357.81 Safari/537.36\n" +
                    "Request Payload\n");
            multipart.addFormField("pdb", row.getPdbCode());
            multipart.addFormField("chain", row.getPdbChain());
            multipart.addFormField("tol", "40");
            multipart.addFormField("len", "20");
            multipart.addFormField("upchainid", "");
            multipart.addFilePart("uploadfile", new File("data/rmsd.o"));
            multipart.addFormField("type1", "1");
            List<String> response = multipart.finish();

//            System.out.println("SERVER REPLIED:");

            for (String line : response) {

                if (line.contains("There is no repeat found for the given parameters")) {
                    //System.out.println(line);
                    //System.out.println(row.getProtein()+"\t"+"0");
                    results = row.getProtein() + "\t" + "0";
                    break;
                } else if (line.contains("second.pl?" + row.getPdbCode() + "&" + row.getPdbChain())) {
                    String url = "http://cluster.physics.iisc.ernet.in/cgi-bin/prostrip/second.pl?" + row.getPdbCode() + "&" + row.getPdbChain();
                    Document doc = Jsoup.connect(url).get();
                    Elements links = doc.select("input[name=options1]");
                    StringBuffer value = new StringBuffer();
                    for (Element link : links) {
                        value.append(link.attr("value") + " ");
                    }

                    //System.out.println(row.getProtein()+"\t"+"1"+"\t["+value.toString()+"]");
                    results = row.getProtein() + "\t" + "1" + "\t[" + value.toString() + "]";
                    break;

                } else {
                    //System.out.println(row.getProtein()+"\t"+"0"+"\t"+"There is no repeat found for the given parameters");
                    results = row.getProtein() + "\t" + "0" + "\t" + "There is no repeat found for the given parameters";
                }
            }
        } catch (IOException ex) {
            System.err.println(ex);
        }
        return results;
        //return urlResult;

    }
}
