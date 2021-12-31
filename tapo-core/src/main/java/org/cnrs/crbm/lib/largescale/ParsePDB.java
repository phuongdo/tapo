package org.cnrs.crbm.lib.largescale;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.CombineScore;

import java.util.List;

/**
 * Created by pdoviet on 6/2/2015.
 */
public class ParsePDB {
    public static void main(String[] args) {

        String fileIn ="data/resultLargeScale/PDBJune2015.o";
        String fileOut = "C:/Users/pdoviet/Desktop/PDB.o";
        new ParsePDB().convertToCSV(fileIn,fileOut);

    }


    public void convertToCSV(String fileIn, String fileOut) {
        List<String> pdbLines = DataIO.readLines(fileIn);
        StringBuffer buffer = new StringBuffer();
        buffer.append("protein,vScore,tmScore,ceScore,caScore,mScore,hScore,sigScore\n");
        for (String line : pdbLines) {
            String pdb = line.split("\t")[0];
            String strScores = line.split("\t")[1];
            CombineScore score = new CombineScore(strScores);
            buffer.append(pdb + "," + score.toCSVFormat() + "\n");

        }

        DataIO.writeToFile(buffer.toString(), fileOut);

    }

}
