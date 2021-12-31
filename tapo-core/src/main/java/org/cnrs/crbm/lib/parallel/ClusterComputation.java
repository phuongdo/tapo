package org.cnrs.crbm.lib.parallel;

import org.cnrs.crbm.lib.conf.ConfigUtil;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.conf.Var;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.utils.ExecUtils;

//import java.io.IOException;
//
//import org.apache.commons.exec.CommandLine;
//import org.apache.commons.exec.DefaultExecutor;
//import org.apache.commons.exec.ExecuteException;

import java.util.Calendar;

/**
 * Created by pdoviet on 4/17/2015.
 * support genotul
 */
public class ClusterComputation {


    String pdbCode = "";
    String pdbChain = "";
    String pid = "";
    int nCore = Var.PARALLEL_SMP;
    String CLUS_MOD = Var.CLUS_MOD;

    String script = "";

    private final static String JAVA_OPT = ConfigUtil.getInstance().getProperty("JAVA_OPT");

    public String getDirOutputFile() {
        return dirOutputFile;
    }

    String dirOutputFile = "";

    public ClusterComputation(String pdbCode, String pdbChain, String pid) {
        this.pdbCode = pdbCode;
        this.pdbChain = pdbChain;
        this.pid = pid;
    }

    public static void main(String[] args) {
        ClusterComputation computation = new ClusterComputation("1lxa", "A", "xxxxxxx");
        computation.generateScript();
//        computation.qsub();
    }


    public void generateScript() {

        //;long jobid = Calendar.getInstance().getTimeInMillis();
        // outscript
        script = Dir.TMP_DIR + "/script_" + pid + "_" + pdbCode + "" + pdbChain + ".sh";
        //Dir.HOME_DIR
        String output = Dir.TMP_DIR + "/" + pid + "_" + pdbCode + "" + pdbChain + ".o";
        //dirOutputFile = output;
        String cmd = "#!/bin/sh\n";
        //cmd += "java -Xmx2042m -Djava.awt.headless=true -cp  \"target/classes;lib/*\" apps.TaPo apps.TaPo -p " + pdbCode + " -c " + pdbChain + " -nCore " + nCore + " -o " + output;
        //String dir = "F:/Code/intelws/tapo";
//        if (nCore >= 5)
//        nCore = nCore / 2;
        cmd += JAVA_OPT + " apps.TaPo -p " + pdbCode + " -c " + pdbChain + " -nCore " + nCore + " -o " + output;
//        if (CLUS_MOD.equals("LOCAL")) {
//            cmd += JAVA_OPT + " apps.TaPo -p " + pdbCode + " -c " + pdbChain + " -nCore " + nCore + " -o " + output;
//        } else
//            cmd += JAVA_OPT + " apps.TaPo -p " + pdbCode + " -c " + pdbChain + " -nCore " + nCore + " -o " + output;
        DataIO.writeToFile(cmd, script);
    }


    public void submitJob() {
        // if local
        String jobName = "job_" + pid + "_" + pdbCode + "" + pdbChain;

        String jobIDOUT = Dir.TMP_DIR + "/" + jobName + ".ID";

        String qsubOption = "-l mem=10G -l h_vmem=15G";
        String shellCommand = "qsub " + qsubOption + " -pe parallel_smp " + Var.PARALLEL_SMP + " -o " + Dir.TMP_DIR + "/" + " -e " + Dir.TMP_DIR + "/" +
                " " + script + " | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}' > " + jobIDOUT;
        //System.out.println(shellCommand);
        //this.runScript(shellCommand);
        //System.out.println(strCmdOut);
        if (CLUS_MOD.equals("LOCAL"))
            shellCommand = script;
        ExecUtils.execShellCmdLinux(shellCommand);

    }


}
