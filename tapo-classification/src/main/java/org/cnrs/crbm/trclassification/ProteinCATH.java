package org.cnrs.crbm.trclassification;

/**
 * Created by pdoviet on 9/16/2015.
 */
public class ProteinCATH {


    public String getCathId() {
        return cathId;
    }

    public void setCathId(String cathId) {
        this.cathId = cathId;
    }

    public String getAssignCathId() {
        return assignCathId;
    }

    public void setAssignCathId(String assignCathId) {
        this.assignCathId = assignCathId;
    }

    //    private String cathId = "";
    private String cathId = "unk";
    private String assignCathId = "4";

    @Override
    public String toString() {
        return assignCathId + "\t" + cathId;
    }
}
