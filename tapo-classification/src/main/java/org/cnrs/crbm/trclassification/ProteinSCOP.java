package org.cnrs.crbm.trclassification;

/**
 * Created by pdoviet on 9/16/2015.
 */
public class ProteinSCOP {

    public String getScopId() {
        return scopId;
    }

    public void setScopId(String scopId) {
        this.scopId = scopId;
    }

    public String getAssignScopId() {
        return assignScopId;
    }

    public void setAssignScopId(String assignScopId) {
        this.assignScopId = assignScopId;
    }

    //    private String cathId = "";
    private String scopId = "unk";
    private String assignScopId = "x";

    @Override
    public String toString() {
        return assignScopId + "\t" + scopId;
    }
}
