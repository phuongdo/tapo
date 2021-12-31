package org.cnrs.crbm.trclassification;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.net.URL;
import java.net.URLConnection;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import com.apporiented.algorithm.clustering.*;
import org.apache.commons.io.FileUtils;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.scop.ScopDomain;
import org.biojava.nbio.structure.scop.ScopInstallation;
import org.cnrs.crbm.lib.classification.JobDistanceCals;
import org.cnrs.crbm.lib.classification.ProData;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.io.*;
import org.cnrs.crbm.lib.math.MapUtil;
import org.cnrs.crbm.lib.parallel.JobInput;
import org.cnrs.crbm.lib.repeats.Superimposer;
import org.cnrs.crbm.lib.repeats.SuperimposerOutput;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

public class ProClassify {

    static final String JDBC_DRIVER = "com.mysql.jdbc.Driver";
    // using ssh tunning port forwarding
    static final String DB_URL = "jdbc:mysql://localhost:3306/pbank";

    // Database credentials
    static final String USER = "pdoviet";
    static final String PASS = "kajava";

    public static void main(String[] args) throws Exception {

        // Protein Domain Parser (PDP)

        // new ProClassify().proteinTRsAnalysis();
        // new ProClassify().classify();
        //new ProClassify().hierarchicalClustering();
        //new ProClassify().saveCluster();
        // new ProClassify().saveCluster();
        // System.out.println("PDP:2IF4Aa".substring(4, "PDP:2IF4Aa".length()));
        // new ProClassify().saveCluster();
        // System.out.println(new ProClassify().getRepresentativeDomains("3v6o",
        // "C"));

        // new ProClassify().threads();

//		new ProClassify().generateFilesForLargeScaleAnalysis();
        new ProClassify().hierarchicalClustering();
//		new ProClassify().classifyProteinByScope("1yi2","G");
//        System.out.println(new ProClassify().classifyProteinBySimilarity("4akg", "A"));

    }


    private void generateFilesForLargeScaleAnalysis()
            throws FileNotFoundException {
        PrintWriter writer = new PrintWriter("data/dataset1.clus.in");
        List<Row> rows = DataIO.getListProteinFromFile("data/benchmarkset/dataset1.in");
        for (int i = 0; i < rows.size(); i++) {
            Row rowi = rows.get(i);
            for (int j = i + 1; j < rows.size(); j++) {
                Row rowj = rows.get(j);

                writer.write(rowi.getProtein() + "\t" + rowj.getProtein()
                        + "\n");

            }
        }

        writer.close();

    }

    private void threads() {
        int nrOfProcessors = 6;
        ExecutorService eservice = Executors.newFixedThreadPool(nrOfProcessors);
        List<Row> rows = DataIO.getListProteinFromFile("data/benchmarkset/dataset1.in");
        JobInput jobInput = new JobInput();
        jobInput.setProteins(rows);

//		int jobid = 0;
        for (int jobid = 0; jobid < rows.size(); jobid++) {

            Row row = rows.get(jobid);
            eservice.execute(new JobDistanceCals(jobid, jobInput, row));
            //jobid++;
        }

        if (eservice.isTerminated()) {

        }

        eservice.shutdown();

    }

    List<RowRepeatDB> rowsRef = null;
    Superimposer superimpose = new Superimposer();
    Map<String, String> pdbTRs = new HashMap<String, String>();
    final Taxonomy t1 = new Taxonomy();
    final Taxonomy t2 = new Taxonomy();

    public ProClassify() throws Exception {

        // setup
        this.rowsRef = this.getRepeatDB();
//        t1.setScientificName("TRs");
//        t2.setScientificName("No-TRs");
//        // read the output from TAPO analyses.(rossmann.TRs.o)
//        ReadFasta fasta = new ReadFasta();
//        Map<String, String> a = fasta.getFastaFileInSpecialCase("data/rossmann/rossmann.TRs.o");
//        for (Map.Entry<String, String> entry : a.entrySet()) {
////            System.out.println(entry.getKey());
//            //System.out.println(">>>>>>" + entry.getValue());
//            String pdb = entry.getKey().split("\\|")[0];
//            String pdbCode = pdb.substring(0, 4);
//            String pdbChain = pdb.substring(5, 6);
//            String scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, 0, 0);
//
//            if (entry.getValue().contains("No-TRs")) {
//                pdbTRs.put(pdb, "No-TRs");
//            } else {
//                pdbTRs.put(pdb, "3DTRs");
//            }
//
//        }
    }

    public void proteinTRsAnalysis() throws Exception {
        List<String> lines = DataIO.readLines("data/NouvelleTRs.txt");
        List<String> domainLines = DataIO.readLines("data/domain.cluster");
        Map<String, String> domMap = new HashMap<String, String>();
        String dirOut = Dir.CLUSTER_WORKING_DIR + "/output";
        System.out.println("deleting all file and folder in  : " + dirOut);
        PrintWriter writer = new PrintWriter("data/pdb.cluster");
        FileUtils.cleanDirectory(new File(dirOut));

        for (String l : domainLines) {

            domMap.put(l.split("\\;")[0], l.split("\\;")[1]);

        }

        String line = null;

        List<ProData> listpros = new ArrayList<ProData>();
        // ignore header
        for (int i = 1; i < lines.size(); i++) {
            line = lines.get(i);

            String[] rows = line.split("\t");
            String protein = rows[0];
            String strClass = rows[1];
            int nrTRs = Integer.parseInt(rows[2]);
            String anno = rows[3];
            int len = Integer.parseInt(rows[4]);
            // System.out.println(rows[5]);

            ProData prodata = new ProData();
            prodata.setAnno(anno);
            prodata.setProtein(protein);
            prodata.setLen(len);
            prodata.setNrTRs(nrTRs);
            prodata.setStrClass(strClass);
            if (!rows[5].equals("[]")) {

                String[] domains = rows[5].substring(1, rows[5].length() - 1)
                        .split("\\,");
                for (String domain : domains) {
                    prodata.getDomains().add(domain.trim());

                }
            }

            listpros.add(prodata);

        }

        Connection conn = null;
        Statement stmt = null;

        // STEP 2: Register JDBC driver
        Class.forName("com.mysql.jdbc.Driver");

        // STEP 3: Open a connection
        System.out.println("Connecting to database...");
        conn = DriverManager.getConnection(DB_URL, USER, PASS);

        stmt = conn.createStatement();

        int i = 1;
        for (ProData data : listpros) {

            StringBuffer buffer = new StringBuffer();
            for (String domain : data.getDomains()) {

                if (domMap.containsKey(domain)) {
                    buffer.append(domMap.get(domain) + "&");
                } else {
                    // compare
                    for (Entry<String, String> entry : domMap.entrySet()) {
                        try {
                            if (isSimilarDomain(domain, entry.getKey(), stmt)) {
                                buffer.append(entry.getValue() + "&");

                            }
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }

                }

            }

            String classification = "null";
            String protein = data.getProtein();

            if (buffer.length() == 0) {
                buffer.append("none\tnone");

                writer.write(data.getProtein() + "\t" + buffer.toString()
                        + "\n");
                classification = "none";

            } else {

                String[] classes = buffer.toString().split("\\&");

                Map<String, Integer> freqs = new HashMap<String, Integer>();

                for (String strClass : classes) {

                    if (freqs.containsKey(strClass)) {
                        freqs.put(strClass, freqs.get(strClass) + 1);
                    } else {
                        freqs.put(strClass, 1);
                    }

                }

                freqs = MapUtil.sortByValue(freqs);
                classification = (String) freqs.keySet().toArray()[0];
                writer.write(data.getProtein() + "\t"
                        + freqs.keySet().toArray()[0] + "\t" + freqs + "\n");
            }

            // new File(dirOut + "/cluster" + i).mkdirs();

            String directory = dirOut + "/" + classification;
            String pdbFile = directory + "/" + protein + ".pdb";
            if (!new File(directory).exists()) {

                new File(directory).mkdir();
            }
            // save proteins;

            String pdbCode = protein.substring(0, 4);
            String pdbChain = protein.substring(5, 6);

            Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);

            Chain c = structure.getChainByPDB(pdbChain);
            Structure newStruct = new StructureImpl();
            newStruct.addChain(c);

            DataIO.writeToFile(newStruct.toPDB(), pdbFile);
            System.out.print("\r" + i + "/" + listpros.size());

            i++;
        }
        writer.close();
        stmt.close();
        conn.close();
        // Map<String, String> map = new HashMap<String, String>();
        //
        // for (ProData data : listpros) {
        //
        // for (String domain : data.getDomains()) {
        //
        // if (map.containsKey(domain)) {
        //
        // map.put(domain, map.get(domain) + "," + data.getProtein());
        //
        // } else {
        // map.put(domain, data.getProtein());
        // }
        //
        // }
        //
        // }
        //
        // for (Entry<String, String> entry : map.entrySet()) {
        //
        // System.out.println(entry.getKey() + "\t" + entry.getValue());
        // }

    }

    private boolean isSimilarDomain(String domain, String otherDomain,
                                    Statement stmt) throws SQLException {

        if (domain.equals(otherDomain))
            return true;
        String sql = this.buildQuery(domain, otherDomain);
        // System.out.println(sql);
        ResultSet rs = stmt.executeQuery(sql);
        while (rs.next()) {
            // Retrieve by column name
            int sim1 = rs.getInt("sim1");
            int sim2 = rs.getInt("sim2");
            // int len1 = rs.getInt("len1");
            // int len2 = rs.getInt("len2");
            int sim_max = Math.max(sim1, sim2);
            int sim_min = Math.min(sim1, sim2);
            if (sim_max > 80 && sim_min > 40)
                return true;
        }
        return false;
    }

    final static String RDOMAIN_URL = "http://www.rcsb.org/pdb/rest/representativeDomains?structureId=";

    public void saveCluster() throws IOException, StructureException {
        List<String> rows = DataIO.readLines("data/pdb.hire.cluster");
        System.out.println("cleanning output...");
        String dirOut = Dir.CLUSTER_WORKING_DIR + "/output";
        FileUtils.cleanDirectory(new File(dirOut));

        for (int i = 0; i < rows.size(); i++) {
            new File(dirOut + "/cluster" + i).mkdirs();

            System.out.println(rows.get(i));
            String[] domainIDs = rows.get(i).split(" ")[1].split("\\&");
            for (String domainid : domainIDs) {

                try {
                    // Structure domain = StructureIO.getStructure(domainid);

                    String pdbCode = domainid.substring(0, 4);
                    String pdbChain = domainid.substring(5, 6);

                    Structure structure = PdbTools
                            .getStructureFromLocalPdb(pdbCode);

                    Chain c = structure.getChainByPDB(pdbChain);
                    Structure domain = new StructureImpl();
                    domain.addChain(c);

                    String fileDir = dirOut + "/cluster" + i + "/" + domainid
                            + ".pdb";

                    System.out.println(fileDir);
                    DataIO.writeToFile(domain.toPDB(), fileDir);

                } catch (Exception ex) {

                    ex.printStackTrace();
                }
            }
        }

        // Structure domain1 = StructureIO.getStructure("PDP:3IQUAb");

        // save structure.

        // DataIO.writeToFile(domain1.toPDB(), "output/3iquab.pdb");
    }

    public List<String> getRepresentativeDomains(String pdbCode, String pdbChain)
            throws Exception {

        List<String> pdpDomains = new ArrayList<String>();

        URL url = new URL(RDOMAIN_URL + pdbCode + "." + pdbChain);

        URLConnection connection = url.openConnection();

        Document doc = parseXML(connection.getInputStream());
        NodeList descNodes = doc.getElementsByTagName("representative");

        for (int i = 0; i < descNodes.getLength(); i++) {
            // System.out.println(descNodes.item(i).getTextContent());
            Element el = (Element) descNodes.item(i);
            String domain = el.getAttribute("name");
            pdpDomains.add(domain);
        }

        return pdpDomains;
    }

    private Document parseXML(InputStream stream) throws Exception {
        DocumentBuilderFactory objDocumentBuilderFactory = null;
        DocumentBuilder objDocumentBuilder = null;
        Document doc = null;
        try {
            objDocumentBuilderFactory = DocumentBuilderFactory.newInstance();
            objDocumentBuilder = objDocumentBuilderFactory.newDocumentBuilder();

            doc = objDocumentBuilder.parse(stream);
        } catch (Exception ex) {
            throw ex;
        }

        return doc;
    }

    public String classifyProteinByScope(String pdbCode, String pdbChain) {

        StringBuffer strBuffer = new StringBuffer();

        String cacheLocation = Dir.SCOPE_LOCAL;
        String pdbId = "1mhj";
        // download SCOP if required and load into memory
        ScopInstallation scop = new ScopInstallation(cacheLocation);
//        scop.setScopVersion("2.05");
        // scop = (ScopInstallation)
        // ScopFactory.getSCOP(ScopFactory.VERSION_1_75);


        List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);

        for (ScopDomain domain : domains) {
            System.out.println(domain);
            System.out.println(domain.getRanges());
            System.out.println(domain.getClassificationId());
            System.out.println(domain.getScopId());
        }
//		ScopNode node = scop.getScopNode(domains.get(0).getSunid());
//
//		while (node != null) {
//
//			System.out.println("This node: sunid:" + node.getSunid());
//			System.out.println(scop.getScopDescriptionBySunid(node.getSunid()));
//			node = scop.getScopNode(node.getParentSunid());
//		}

        return strBuffer.toString();
    }

    public String classifyProteinBySimilarity(String pdbCode, String pdbChain)
            throws StructureException {

        StringBuffer strBuffer = new StringBuffer();

        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);

        Atom[] atomSet1 = PdbTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain));

        // System.exit(1);
        String assignedBestClass = "UNK";
        double maxTMScore = -1;
        String assignedPdb = "UNK";

        for (RowRepeatDB rowRef : this.rowsRef) {

            try {
                // System.out.print("\Rlib Scan " + rowRef.getEntry() + "/"
                // + rowsRef.size());
                // System.out.println(rowRef.getStrclass());

                Structure structureRef = PdbTools
                        .getStructureFromLocalPdb(rowRef.getPdbCode());

                Atom[] atomSetRef = PdbTools.getAtomCAArray(structureRef
                        .getChainByPDB(rowRef.getPdbChain()));

                // double rmsd = 0.0;
                int len1 = atomSet1.length;
                int len2 = atomSetRef.length;
                // double sim1 = 0.0;

                int avg = (len1 + len2) / 2;

                SuperimposerOutput output = superimpose.compareTwoStructures(atomSet1, atomSetRef);
                if (output.getTmScore() > maxTMScore) {
                    maxTMScore = output.getTmScore();
                    assignedBestClass = rowRef.getStrclass();
                    assignedPdb = rowRef.getPdbCode() + "_" + rowRef.getPdbChain();
                    //strBuffer.append(rowRef.getStrclass());
                }
//                if ((double) Math.abs(avg - len1) / len1 < 0.1)
//
//                {
//
//                    if (superimpose.isTheSameStructure(atomSet1, atomSetRef)) {
//                        // System.out.println(">>" + rowRef.getEntry() +
//                        // " class:"
//                        // + rowRef.getStrclass());
//                        strBuffer.append(rowRef.getStrclass());
//                        break;
//                    }
//
//                }

            } catch (Exception e) {
                // TODO: handle exception
            }
        }

        strBuffer.append(assignedBestClass + "\t" + assignedPdb + "\t" + NumberFormatUtils.format(maxTMScore));

//        if (maxTMScore >= 0.5)
//            strBuffer.append(assignedBestClass);

//        if (strBuffer.length() == 0)
//            strBuffer.append("UNK");
        return strBuffer.toString();

    }

    private void hierarchicalClustering() throws FileNotFoundException {
        List<String> rowsOld = DataIO.readLines("data/rossmann/rossmannPdbList.in");

        List<String> rows = new ArrayList<String>();
        for (String row : rowsOld) {
            //System.out.println(row.split("\t")[0]);
            //rows.add(row.split("\t")[0].trim());
            String pdb = row.split("\t")[0].trim();
            /**
             * for scope trclassification
             */
            String pdbCode = pdb.substring(0, 4);
            String pdbChain = pdb.substring(5, 6);
            ProteinSCOP scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, 0, 0);
//            System.out.println(scopeDomain);
//            if (scopeDomain.equals("c")) {
//                scopeDomain = "a/b";
//            } else if (scopeDomain.equals("d")) {
//                scopeDomain = "a+b";
//            }


            String classTRs = "No-TRs";
            if (pdbTRs.containsKey(pdb)) {
                classTRs = pdbTRs.get(pdb);
            }

            if (classTRs.equals("No-TRs"))
                classTRs = "x";
            else classTRs = "o";

            rows.add(pdb + "_" + scopeDomain + "_" + classTRs);
//            rows.add(pdb);
        }

        //System.exit(1);
        String[] names = this.toName(rows);
        double[][] distances = this.getDistances(rows);

        this.printDistanceMatrix(rows, distances);
        //System.exit(1);
        ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
        Cluster cluster = alg.performClustering(distances, names,
                new CompleteLinkageStrategy());
        showPhylogeneticTree(cluster);
//		traverse(cluster);

        // DendrogramPanel dp = new DendrogramPanel();
        // dp.setModel(cluster);
        //
        // // 1. Create the frame.
        // JFrame frame = new JFrame("FrameDemo");
        // frame.setSize(300, 400);
        //
        // // 2. Optional: What happens when the frame closes?
        // frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        //
        // // 3. Create components and put them in the frame.
        // // ...create emptyLabel...
        // frame.getContentPane().add(dp, BorderLayout.CENTER);
        //
        // // 4. Size the frame.
        // frame.pack();
        //
        // // 5. Show it.
        // frame.setVisible(true);

    }

    private void printDistanceMatrix(List<String> rows, double[][] distances)
            throws FileNotFoundException {

        PrintWriter writer = new PrintWriter("output/distance.o");

        int size = rows.size();
        //writer.write(size + "\n");
        for (int i = 0; i < size; i++) {
            //writer.write(rows.get(i) + "\t");

            for (int j = i + 1; j < size; j++) {

                writer.write(distances[i][j] + "\n");

            }
            //writer.write("\n");

        }

        writer.close();

    }

    private void showPhylogeneticTree(Cluster cluster) {

        Phylogeny phy = new Phylogeny();
        PhylogenyNode root = new PhylogenyNode();
        // PhylogenyNode d1 = new PhylogenyNode();
        // PhylogenyNode d2 = new PhylogenyNode();
        // root.setName("root");
        // d1.setName("descendant 1");
        // d2.setName("descendant 2");
        // root.addAsChild(d1);
        // root.addAsChild(d2);


        buildTree(cluster, root);

        // System.out.println("done");
        phy.setRoot(root);
        phy.setRooted(true);
        // Displaying the newly created tree with Archaeopteryx.
        Archaeopteryx.createApplication(phy);

    }

    public void buildTree(Cluster cluster, PhylogenyNode node) {

        for (Cluster child : cluster.getChildren()) {
            PhylogenyNode childNode = new PhylogenyNode();

            if (!child.isLeaf()) {
                childNode.setName(child.getDistance() + "");
                node.addAsChild(childNode);
                childNode.setDistanceToParent(child.getDistance());
                buildTree(child, childNode);
            } else {
                childNode.setName(child.getName());
//                String pdb = child.getName().substring(0,6);
//                String classTRs = "No-TRs";
//                if(pdbTRs.containsKey(pdb)){
//                    classTRs = pdbTRs.get(pdb);
//                }
//                if(classTRs.equals("No-TRs"))
//                    childNode.getNodeData().addTaxonomy(t2);
//                else
//                    childNode.getNodeData().addTaxonomy(t1);
                node.addAsChild(childNode);
            }

        }

    }

    public void traverse(Cluster cluster) {

        // Each child of a tree is a root of
        // its subtree.
        for (Cluster child : cluster.getChildren()) {
            // if(child.getDistance())
            // if (child.getDistance() < 30) {
            double distance = child.getDistance() == null ? 100.0D : child
                    .getDistance().doubleValue();
            if (distance < 0.5)
                System.out.println(child);
            else
                traverse(child);
            // }

        }

    }

    private String buildQuery(String row1, String row2) {
        String sql;

        String name1 = row1.split("\t")[0];
        String name2 = row2.split("\t")[0];
        // non
        sql = "SELECT `score`,`len1`,`len2`,`rmsd`,`sim1`,`sim2`  FROM `fastcat_align_pdb` WHERE (id='"
                + name1 + "' AND pair = '" + name2 + "')";

        return sql;
    }

    private String buildQuery(Row row1, RowRepeatDB row2) {
        String sql;

        String name1 = row1.getPdbCode().toUpperCase() + row1.getPdbChain();
        String name2 = row2.getPdbCode().toUpperCase() + row2.getPdbChain();
        // non
        sql = "SELECT `score`,`len1`,`len2`,`rmsd`,`sim1`,`sim2`  FROM `fastcat_align_pdb` WHERE (id='"
                + name1
                + "' AND pair = '"
                + name2
                + "') OR (id='"
                + name2
                + "' AND pair = '" + name1 + "')";

        return sql;
    }

    private Map<String, Double> getCompareCache() {

        Map<String, Double> map = new HashMap<String, Double>();

        List<String> lines = DataIO.readLines("data/rossmann/rossmannDistMatrix.cache");

        for (String line : lines) {

            String[] rows = line.split("\t");

            map.put(rows[0] + "_" + rows[1], Double.parseDouble(rows[2]));

        }

        return map;

    }

    private double[][] getDistances(List<String> rows)
            throws FileNotFoundException {

        Connection conn = null;
        Statement stmt = null;

        Map<String, Double> cache = this.getCompareCache();
        double[][] distances = new double[rows.size()][rows.size()];
        try {
            for (int i = 0; i < rows.size(); i++) {

                // String pdbCode = rows.get(i).substring(0, 4);
                // String pdbChain = rows.get(i).substring(5, 6);

                for (int j = 0; j < rows.size(); j++) {

                    double score = 1;
                    if (i == j) {
                        score = 0;

                    } else {
                        try {
                            // pdbCode = rows.get(j).substring(0, 4);
                            // pdbChain = rows.get(j).substring(5, 6);
                            String key1 = rows.get(i).substring(0, 6) + "_" + rows.get(j).substring(0, 6);
                            String keyInverse = rows.get(j).substring(0, 6) + "_" + rows.get(i).substring(0, 6);
//							String key1 = rows.get(i) + "_" + rows.get(j);
//							String keyInverse = rows.get(j) + "_" + rows.get(i);
                            if (cache.containsKey(key1)) {
                                score = cache.get(key1);
                            } else if (cache.containsKey(keyInverse)) {
                                score = cache.get(keyInverse);
                            }

                        } catch (Exception ex) {
                            // System.out.println(pdbCode);
                            ex.printStackTrace();
                        }

                    }

                    // System.out.println(score);
//					if (score < 0.5)
//						score = 0.0;
//					else
//						score = 1;
                    distances[i][j] = score;

                }

            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }
        // writer.close();
        return distances;
    }

    private String[] toName(List<String> rows) {

        String[] names = new String[rows.size()];

        for (int i = 0; i < rows.size(); i++) {
            names[i] = rows.get(i).split("\t")[0].trim();
        }

        return names;

    }

    public void classify() throws Exception {

        // read data from repeatDB

        List<RowRepeatDB> rowsRef = this.getRepeatDB();

        Map<String, String> map = new HashMap<String, String>();

        int i = 0;
        for (RowRepeatDB row : rowsRef) {
            try {
                System.out.print("\r" + i + "/" + rowsRef.size());

                String strClass = row.getStrclass();
                List<String> domains = this.getRepresentativeDomains(
                        row.getPdbCode(), row.getPdbChain());
                for (String domain : domains) {

                    map.put(domain, strClass);
                }

            } catch (Exception ex) {
                ex.printStackTrace();
            }

            i++;
        }

        PrintWriter writer = new PrintWriter("data/domain.cluster");

        //
        for (Map.Entry<String, String> entry : map.entrySet()) {

            // System.out.println(entry.getKey() + ";" + entry.getValue());
            writer.write(entry.getKey() + ";" + entry.getValue() + "\n");
        }
        writer.close();
    }

    List<RowRepeatDB> getRepeatDB() throws Exception {
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");
        List<RowRepeatDB> refs = new ArrayList<RowRepeatDB>();
        for (RowRepeatDB row : rows) {
            if (row.getAnnlevel().equals("Detailed")) {
                refs.add(row);
            }
        }

        return refs;
        // return refs;

    }

}
