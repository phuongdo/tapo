package org.cnrs.crbm.lib.parallel;

import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.conf.Var;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.utils.ExecUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 4/20/2015.
 */
public class ClusterJobStatus {

    String pdbCode = "";
    String pdbChain = "";
    String pid = "";

    public static void main(String[] args) {
        ClusterJobStatus clusterJobStatus = new ClusterJobStatus(args[0], args[1], args[2]);
        System.out.println(clusterJobStatus.getStatus());

    }

    public ClusterJobStatus(String pdbCode, String pdbChain, String pid) {
        this.pdbCode = pdbCode;
        this.pdbChain = pdbChain;
        this.pid = pid;
    }

    public String getStatus() {

        String jobName = "job_" + pid + "_" + pdbCode + "" + pdbChain;
        // get jobID

        String jobIDOUT = Dir.TMP_DIR + "/" + jobName + ".ID";
        String jobID = DataIO.readFile(jobIDOUT).trim();
        String shellCommand = "qstat -u " + Var.USER_NAME + " | grep " + jobID;
        String rep = ExecUtils.execShellCmdLinux(shellCommand);
        //System.out.println(shellCommand);
        //System.out.println(rep);
        // process
        return parseStatus(rep);
    }


    public String parseStatus(String rep) {
        String status = "e";

        if (rep.length() > 5) {
            String[] datas = rep.split(" ");
            List<String> cols = new ArrayList<String>();
            for (String data : datas) {
                //System.out.println(data);
                if (!data.equals(""))
                    cols.add(data);
            }
            status = cols.get(4);
        } else {
            status = "f";
        }
        return status;

    }


}
