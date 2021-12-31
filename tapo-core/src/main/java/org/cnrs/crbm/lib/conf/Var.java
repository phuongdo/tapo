package org.cnrs.crbm.lib.conf;

public class Var {


    public final static int PARALLEL_SMP = Integer.parseInt(ConfigUtil.getInstance()
            .getProperty("PARALLEL_SMP"));

    public final static String CLUS_MOD = ConfigUtil.getInstance()
            .getProperty("CLUSTER_MODE");
    public final static String USER_NAME = ConfigUtil.getInstance()
            .getProperty("USER_NAME");
}
