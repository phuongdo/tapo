package org.cnrs.crbm.nextgen;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Repeat;

import java.util.Map;

/**
 * Program: TAPO v1.1.2: TAPO: A combined method for the identification of tandem repeats in protein structures
 * Authors: Do Viet P, Roche DB, Kajava AV.
 * Report_file: targetList.taf
 * COLUMNS         DATA TYPE           CONTENTS
 * ----------------------------------------------
 * $1              String              PDB_ID
 * $2              String              Finder Name
 * $3              String              Cluster name(selected mark is the best prediction)
 * $4              Integer             Number of repeat units
 * $5              Real                Avg length of repeat
 * $6              String              TRs regions
 * $7              String              Repeat Units regions
 * $8              Real                QA score
 * $9              Real                R Score
 */


public class ConvertTaPoFormat {

    public static void main(String[] args) throws Exception {
        String fileDir = "data/francois/italy.out";
        String fileDirConverterd = "data/francois/italy.out.converted";
        ConvertTaPoFormat convertTaPoFormat = new ConvertTaPoFormat();
        convertTaPoFormat.convert(fileDir, fileDirConverterd);
    }

    public void convert(String fileDir, String fileDirConverterd) throws Exception {

        StringBuffer buffer = new StringBuffer();
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat(fileDir);
        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
            TaPoFastaFormat taPoFastaFormat = entry.getValue();
            String pdbId = taPoFastaFormat.getPdbID();
            String pdbCode = pdbId.substring(0, 4);
            String pdbChain = pdbId.substring(5, 6);
            //extract atoms
            RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
            repeatFinder.setOutputRelativeScore(false);
            buffer.append(">" + taPoFastaFormat.getHeader() + "\n");
            if (taPoFastaFormat.is3DRepeat()) {
                for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                    buffer.append(this.convertToPymolVisualizationFormat(repeatFinder, repeat) + "\n");
                }
            } else {
                buffer.append("No-TRs\n");
            }
        }

        // write to file

        DataIO.writeToFile(buffer.toString(), fileDirConverterd);

    }


    public String convertToPymolVisualizationFormat(RepeatFinder repeatFinder, Repeat repeat) {
        return repeatFinder.buildOutputOneRepeat(repeat).toString();
    }


}
