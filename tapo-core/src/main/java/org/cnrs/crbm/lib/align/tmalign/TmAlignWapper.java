package org.cnrs.crbm.lib.align.tmalign;

/**
 * Created by pdoviet on 25/03/2015 and  modified in 7/30/2015.
 * A simple wrapper for TM-pairAlign by Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710
 * unpublished, just use for local analysis. His code is modified base on java source of Daniel ROCHE
 */

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.utils.ExecUtils;

import java.io.*;
import java.util.*;
import java.lang.*;

public class TmAlignWapper {
    //    Runtime r = Runtime.getRuntime();
//    Process p = null;
    float TMscore1 = 0;
    float TMscore2 = 0;
    StringBuffer alingbuf = new StringBuffer();
    int alignedLength = 0;
    int targetLength = 0;
    int templateLength = 0;
    float RMSD = 0;

    public TmAlignWapper(String pdbfile1, String pdbfile2, String supfilename) {
        try {

            String bin_path = Dir.TMALIGN_EXECUTABLE;
//            String bin_path = "$TMALIGN_HOME/TMalign";
            String supoutstr = "";
            if (!supfilename.equals("NO_FILE"))
                supoutstr = " -o " + supfilename;
//            String[] runtmalign = {"/bin/bash", "-c", "ulimit -t 3600; cd " + direct + " ;" + bin_path + " " + pdbfile1 + " " + pdbfile2 + supoutstr};
//            p = r.exec(runtmalign);
            String cmd = bin_path + " -A " + pdbfile1 + " -B " + pdbfile2 + supoutstr;
            String[] lines = ExecUtils.execShellCmdLinux(cmd).split("\n");


            boolean appendalingment = false;
            for (String line : lines) {
//                System.out.println(line);
                if (line.startsWith("Length of structure A:")) {
                    String targetlengthstr = line.substring(line.indexOf(":") + 1, line.lastIndexOf("residues"));
                    targetLength = (new Integer(targetlengthstr.trim())).intValue();
                }

                if (line.startsWith("Length of structure B")) {
                    String templatelengthstr = line.substring(line.indexOf(":") + 1, line.lastIndexOf("residues"));
                    templateLength = (new Integer(templatelengthstr.trim())).intValue();
                }


                if (line.startsWith("Aligned length=")) {
                    StringTokenizer linetoks = new StringTokenizer(line, ",");

                    String alignedlengthstr = linetoks.nextToken();
                    alignedLength = (new Integer(alignedlengthstr.substring(alignedlengthstr.indexOf("=") + 1).trim())).intValue();

                    String RMSDstr = linetoks.nextToken();
                    RMSD = (new Float(RMSDstr.substring(RMSDstr.indexOf("=") + 1).trim())).floatValue();
                }

                if (line.startsWith("TM-score=")) {
                    StringTokenizer linetoks = new StringTokenizer(line);

                    if (line.contains("(if normalized by length of structure A")) {
//                        linetoks.nextToken();
//                        String TMscorestr = linetoks.nextToken();
//                        TMscore1 = (new Float(TMscorestr.trim())).floatValue();
                        TMscore1 = (new Float(line.substring(9, 15))).floatValue();

                    }
                    if (line.contains("(if normalized by length of structure B")) {
//                        linetoks.nextToken();
//                        String TMscorestr = linetoks.nextToken();
                        TMscore2 = (new Float(line.substring(9, 15))).floatValue();
                    }
                }

                if (line.startsWith("(\":\" denotes residue pairs of d <  5.0 Angstrom, \".\" denotes other aligned residues)"))
                    appendalingment = true;

                if (appendalingment)
                    alingbuf.append(line + "\n");

            }
        } catch (
                Exception e
                )

        {
            System.out.println("Error executing TMalign!" + e);
        }

    }

    public float getTMscore1() {
        return TMscore1;
    }

    public float getTMscore2() {
        return TMscore2;
    }

    public float getRMSD() {
        return RMSD;
    }

    public int getAlignedLength() {
        return alignedLength;
    }

    public int getTargetLength() {
        return targetLength;
    }

    public int getTemplateLength() {
        return templateLength;
    }

    public String getAlignment() {
        return alingbuf.toString();
    }

    public static void main(String args[]) {
        String top_model = "F:/Cygwin/home/pdoviet/BioApps/TMalignc/data/obj01.pdb";
        String templatePDB = "F:/Cygwin/home/pdoviet/BioApps/TMalignc/data/obj02.pdb";

        TmAlignWapper rtm3 = new TmAlignWapper(top_model, templatePDB, "NO_FILE");
        System.out.println("Target length: " + rtm3.getTargetLength() + "\nTemplate length: " + rtm3.getTemplateLength() + "\nAligned length: " + rtm3.getAlignedLength() + "\nRMSD: " + rtm3.getRMSD() + "\nTM-score target: " + rtm3.getTMscore1() + "\nTM-score template: " + rtm3.getTMscore2() + "\nAlignment:\n" + rtm3.getAlignment());
    }
}
