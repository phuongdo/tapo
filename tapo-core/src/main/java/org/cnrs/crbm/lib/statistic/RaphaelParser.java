package org.cnrs.crbm.lib.statistic;

import org.cnrs.crbm.lib.io.DataIO;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;

import java.io.File;
//import java.util.DoubleSummaryStatistics;
import java.util.List;

/**
 * Created by pdoviet on 5/11/2015.
 */
public class RaphaelParser {

    public static void main(String[] args) throws Exception {

        new RaphaelParser().parse();

    }


    public void parse() {

        List<String> rows = DataIO.readLines("data/tapo/trainsetPdb.in");

        for (String row : rows) {

            String pdbCode = row.substring(0, 4);
            String pdbChain = row.substring(5, 6);
            String strClass = row.split("\t")[1];
            if (strClass.equals("1"))
                strClass = "1";
            else
                strClass = "0";

            if(strClass.equals("1"))
                continue;
            String entry = pdbCode + "_" + pdbChain;

            // get smv score

            double svmScore = -1.0;


            // check html file
            String filePathString = "output/raphael_html/" + entry + ".html";
            File f = new File(filePathString);
            if (f.exists() && !f.isDirectory()) {
            /* extract smv score */
                String html = DataIO.readFile(filePathString);
                Document doc = Jsoup.parse(html);

                //*[@id="inputdata"]/table[1]/tbody/tr/td[2]/pre/font[3]
                Elements svmE = doc.select("#inputdata > table:eq(1) > tbody > tr > td:eq(1) > pre > font:eq(2) ");
                //System.out.println(svmE.text());
                try {
                    String svmStr = svmE.text().substring(11, 20);//.replaceAll(".","").trim();
                    //svmStr = svmStr.replaceAll(".","");
                    //System.out.println(entry + "\t" + svmStr);
                    svmScore = Double.parseDouble(svmStr);
                }catch (Exception ex){

                }
                System.out.println(entry + "\t" + svmScore);

            }


        }


    }


}
