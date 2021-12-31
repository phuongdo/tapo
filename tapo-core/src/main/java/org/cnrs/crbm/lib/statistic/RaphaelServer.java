package org.cnrs.crbm.lib.statistic;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.networking.MultipartUtility;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Created by pdoviet on 2/5/2015.
 */
public class RaphaelServer {
    private static String WEBSERVER = "http://protein.bio.unipd.it/raphael/Raphael.jsp";

    public static void main(String[] args) throws Exception {

        RaphaelServer raphael = new RaphaelServer();
        //raphael.runEvaluation();
        raphael.export();
//
    }


    private void export() {

        List<String> lines = DataIO.readLines("output/raphael/job.tab");

        Set<String> set = new HashSet<String>();
        for (String line : lines) {

            if (line.length() > 0)
                set.add(line);

        }

        for (String line : set) {
            String pdb = line.split("\t")[0];
            String links = line.split("\t")[1];

            downloadResultFromServer(links, pdb);
        }


    }

    public void downloadResultFromServer(String link, String fileName) {

        String fileDir = "output/raphael_html/" + fileName + ".html";
        // check file exist
        if (!new File(fileDir).exists()) {
            // load
            try {

                URL url = new URL(link);
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

    private String submitQuery(Row row) throws Exception {

        String charset = "UTF-8";
        String requestURL = "http://protein.bio.unipd.it/raphael/Raphael.jsp";

        String urlResult = "";
        try {
            MultipartUtility multipart = new MultipartUtility(requestURL, charset);

            multipart.addHeaderField("User-Agent", "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/40.0.2214.94 Safari/537.36");

            multipart.addFormField("pdbcode", row.getPdbCode());
            multipart.addFormField("chain", row.getPdbChain());
            multipart.addFilePart("pdbfile", new File("data/rmsd.o"));

            List<String> response = multipart.finish();

            // System.out.println("SERVER REPLIED:");

            for (String line : response)
                if (line.contains("URL=http://protein.bio.unipd.it/raphael/work")) {
                    //System.out.println(line);
                    Document doc = Jsoup.parse("<html><head>" + line + "</head></html>");
                    String content = doc.getElementsByTag("meta").get(0).attr("content");
                    urlResult = content.split("=")[1];
                    break;
                }
        } catch (IOException ex) {
            System.err.println(ex);
        }
        return urlResult;

    }

    private void runEvaluation() throws Exception {


        Random random = new Random();
        File folder = new File("output/raphael/");
        File[] listOfFiles = folder.listFiles();
        Set<String> set = new HashSet<String>();
        for (File file : listOfFiles) {
            if (file.isFile()) {
                set.add(file.getName());
            }
        }
        List<Row> rows = DataIO.getListProteinFromFile("data/pdbList.in");

        for (Row row : rows) {
            StringBuffer buffer = new StringBuffer();
            // read pdbList

            //System.out.println(row.getProtein());
            String filePathString = "output/raphael/" + row.getProtein() + ".tab";
            File f = new File(filePathString);
            if (set.contains(row.getProtein() + ".tab"))
                continue;

            String urlResults = this.submitQuery(row);
            buffer.append(row.getProtein() + "\t" + urlResults + "\n");
            System.out.println(row.getProtein() + "\t" + urlResults);

            if (urlResults.equals(""))
                break;

            if (f.exists() && !f.isDirectory()) {
                continue;
            } else
                DataIO.writeToFile(buffer.toString(), "output/raphael/" + row.getProtein() + ".tab");

            long minutes = 3 * 60 * 1000;
            long randomSecond = random.nextInt(3 * 60) * 1000;
            Thread.sleep(minutes + randomSecond);
            //break;
        }


    }

}
