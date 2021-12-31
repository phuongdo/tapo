package org.cnrs.crbm.lib.services;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.cdhit.CDHIT;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.conf.Var;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.HtmlGenerator;
import org.cnrs.crbm.lib.parallel.ClusterComputation;
import org.cnrs.crbm.lib.parallel.ClusterJobStatus;
import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

public class HandlerResource extends AbstractHandler {

    static Logger logger = LoggerFactory.getLogger(HandlerResource.class);

    public void handle(String target, Request baseRequest,
                       HttpServletRequest request, HttpServletResponse response)
            throws IOException, ServletException {
        response.setContentType("text/html;charset=utf-8");

        try {
            response.setStatus(HttpServletResponse.SC_OK);
            baseRequest.setHandled(true);

            String queryString = request.getQueryString();
            Map<String, String> parameters = WebUtils
                    .getParameters(queryString);
            // System.out.println(queryString);

            //logger.info(queryString);

            // Now you can get your parameter:
            String service = parameters.get("service");

            if (service.equals("read")) {
                String code = parameters.get("code");
                String chain = parameters.get("chain");
                String pid = parameters.get("pid");
                response.getWriter().println(this.checkResults(code, chain, pid));

//                response.getWriter().println(
//                        "PDB CODE : <b>" + code + "</b> Chain:<b>" + chain
//                                + "</b><br/>");
//                String conString = DataIO.readFile("output/"
//                        + code.toLowerCase() + chain + ".out");
//                response.getWriter().println(conString);

            } else if (service.equals("show")) {
                String code = parameters.get("code");
                String chain = parameters.get("chain");
                String pid = parameters.get("pid");
//                String conString = DataIO.readFile("output/"
//                        + code.toLowerCase() + chain + ".out.combined");
                response.getWriter().println(this.showResutls(code, chain, pid));
//                response.getWriter().println(
//                        "PDB CODE : <b>" + code + "</b> Chain:<b>" + chain
//                                + "</b><br/>");
//
//                response.getWriter().println(conString);


            } else if (service.equals("showAll")) {
                String code = parameters.get("code");
                String chain = parameters.get("chain");
                String pid = parameters.get("pid");

                response.getWriter().println(this.showResutlsAll(code, chain, pid));


            } else if (service.equals("analyse")) {

                String code = parameters.get("code").toLowerCase();
                String chain = parameters.get("chain").toUpperCase();
                String pid = parameters.get("pid");
//                System.out.println("run:"+ code);

                response.getWriter().println(this.submitJob(code, chain, pid));
                //response.getWriter().println(this.callProgram(code));
            }
        } catch (
                Exception ex
                )

        {
            response.getWriter().println("There is something wrong!!!");
        }

    }

    private String showResutls(String pdbCode, String pdbChain, String pid) {

        HtmlGenerator generator = new HtmlGenerator(pdbCode, pdbChain, pid);
        generator.setSHOW_OPT("TheBest");
        generator.buildHTML();
        return generator.html();

    }

    private String showResutlsAll(String pdbCode, String pdbChain, String pid) {

        HtmlGenerator generator = new HtmlGenerator(pdbCode, pdbChain, pid);
        generator.buildHTML();
        return generator.html();

    }

    private String submitJob(String pdbCode, String pdbChain, String pid) throws StructureException {
        ClusterComputation computation = new ClusterComputation(pdbCode, pdbChain, pid);
        computation.generateScript();
        computation.submitJob();
        return "1";

    }

    private String checkResults(String pdbCode, String pdbChain, String pid) {
        String status = "e";
        if (!Var.CLUS_MOD.equals("LOCAL")) {
            ClusterJobStatus clusterJobStatus = new ClusterJobStatus(pdbCode, pdbChain, pid);
            status = clusterJobStatus.getStatus();
        } else {
            status = "f";
        }
        if (status.equals("f")) {
            String dirOutput = Dir.TMP_DIR + "/" + pid + "_" + pdbCode + "" + pdbChain + ".o";
            File f = new File(dirOutput);
            if (!f.exists())
                return "e";
            else {
                String noTRs = pdbCode + "_" + pdbChain + "\tNo-TRs";
                String content = DataIO.readFile(dirOutput).trim();
                if (content.contains(noTRs))
                    return "f_N";
                else return "f_Y";
            }

        }
        return status;


    }

    private String callProgram(String pdbCode) {

        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        CDHIT cdhit = new CDHIT(structure, 0.9, 5);
        List<String> chains = cdhit.getUniqueChain();
        StringBuilder builder = new StringBuilder();
        boolean isRepeat = false;
        boolean tooLargeProtein = false;
        if (chains.size() > 0) {
            // System.out.println("Analyse protein : " + pdbCode);
            for (String pdbChain : chains) {
                try {
                    RepeatFinder repeatFinder = new RepeatFinder(pdbCode,
                            pdbChain);
                    //System.out.println(pdbChain);
                    if (repeatFinder.getAtoms().length > 1000) {
                        tooLargeProtein = true;
                        throw new Exception(
                                "[WARNING] "
                                        + pdbCode
                                        + "_"
                                        + pdbChain
                                        + "  Sorry,  please submit the structure with number of atoms smaller than 1000 !!!");

                    }
                    repeatFinder.setMode(FinderMode.HTML);
                    repeatFinder.findRepeat();
                    if (repeatFinder.isRepeat()) {
                        builder.append(repeatFinder.getOutput() + "<br/>");
                        isRepeat = true;
                    }
                    //System.out.println("DONE");
                    //System.out.println(builder.toString());
                } catch (Exception ex) {

                    builder.append(ex.getMessage() + "</br>");

                    // ex.printStackTrace();
                }
            }
            if (!isRepeat && !tooLargeProtein)
                builder.append("Sorry! No repeat was found!!!");

        } else
            builder.append("CD-HIT throws an exception!!!");
        return builder.toString();

    }
}