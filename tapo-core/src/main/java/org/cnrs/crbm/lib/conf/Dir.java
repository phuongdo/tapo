package org.cnrs.crbm.lib.conf;

public class Dir {

    /**
     * HOME
     */
    public static final String HOME_DIR = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("HOME_DIR");
    public static final String TMP_DIR = ConfigUtil.getInstance().getProperty(
            "TMP_DIR");

    /**
     * CONTACT MAP
     */

    public static final String CM_OUTPUT_DIR = ConfigUtil.getInstance()
            .getProperty("CM_OUTPUT_DIR");

    /**
     * DALILITE 3.3
     */

    public static final String DALI_DIR = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("DALI_DIR");
    public static final String DALI_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("DALI_TMP_DIR");
    /**
     * DSSP
     */
    public static final String DSSP_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("DSSP_EXECUTABLE");
    public static final String DSSP_EXECUTABLE_v2 = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("DSSP_EXECUTABLE_v2");

    public static final String DSSP_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("DSSP_TMP_DIR");

    /**
     * store pdb or fasta file after download
     */
    public final static String DB_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("DB_TMP_DIR");

    /**
     * Structure alphabet temporary directory
     */
    public final static String SADB_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("SADB_TMP_DIR");
    public final static String SADB_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("SADB_EXECUTABLE");

    /**
     * CD-HIT application
     */
    public final static String CDHIT_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("CDHIT_TMP_DIR");
    public final static String CDHIT_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("CDHIT_EXECUTABLE");

    /**
     * MUSTANG (v3.2.1): A MUltiple STuructural AligNment alGorithm.
     */
    public final static String MUSTANG_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("MUSTANG_TMP_DIR");
    public final static String MUSTANG_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("MUSTANG_EXECUTABLE");


    /**
     * 3DCOMB Multiple Protein Structure Alignmen
     */
    public final static String DCOMB_TMP_DIR = ConfigUtil.getInstance()
            .getProperty("DCOMB_TMP_DIR");
    public final static String DCOMB_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("DCOMB_EXECUTABLE");


    /**
     * TM-Align structure alignment program (Zhang and Skolnick 2003).
     * http://zhanglab.ccmb.med.umich.edu/TM-align/
     */

    public final static String TMALIGN_EXECUTABLE = Path.CYGWIN_PATH
            + ConfigUtil.getInstance().getProperty("TMALIGN_EXECUTABLE");

    /**
     * JMOL WEB SITE
     */

    public final static String JMOL_WEB_DATA = ConfigUtil.getInstance()
            .getProperty("JMOL_WEB_DATA");

    /**
     * LARGE-SCALE ANALYSIS
     */
    public final static String PDB_LOCAL = ConfigUtil.getInstance()
            .getProperty("PDB_LOCAL");

    public final static String FASTA_LOCAL = ConfigUtil.getInstance()
            .getProperty("FASTA_LOCAL");

    public final static String DSSP_LOCAL = ConfigUtil.getInstance()
            .getProperty("DSSP_LOCAL");

    public final static String SADB_LOCAL = ConfigUtil.getInstance()
            .getProperty("SADB_LOCAL");

    public final static String SCOPE_LOCAL = ConfigUtil.getInstance()
            .getProperty("SCOPE_LOCAL");
    public final static String CATH_LOCAL = ConfigUtil.getInstance()
            .getProperty("CATH_LOCAL");

    public final static String CLUSTER_WORKING_DIR = ConfigUtil.getInstance()
            .getProperty("CLUSTER_WORKING_DIR");

}
