package org.cnrs.crbm.lib.conf;

/**
 * Created by pdoviet on 3/2/2015.
 */
public class ThresholdConfig {

    public final static double TM_THRES = Double.parseDouble(ConfigUtil.getInstance()
            .getProperty("TM_THRES"));
    public final static double TM_THRES_EXTEND = Double.parseDouble(ConfigUtil.getInstance()
            .getProperty("TM_THRES_EXTEND"));
    public final static double QA_THRES = Double.parseDouble(ConfigUtil.getInstance()
            .getProperty("QA_THRES"));
    public final static double VECTOR_THRES = Double.parseDouble(ConfigUtil.getInstance()
            .getProperty("VECTOR_THRES"));


}
