package org.cnrs.crbm.lib.services;

import org.cnrs.crbm.lib.conf.ConfigUtil;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.ContextHandler;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.handler.ResourceHandler;

/**
 * http://localhost:8888/?code=1PDB&chain=A
 *
 * @author phuongdo
 */
public class WebServer {


    public static void main(String[] args) {
//        org.eclipse.jetty.util.log.Log.setLog(new ServerLogger());
        // org.eclipse.jetty.util.log.Log.setLeve

        int open_port = Integer.parseInt(ConfigUtil.getInstance()
                .getProperty("WEB_PORT"));
        Server server = new Server(open_port);
        try {

            System.out.println("server is ready!!");
            // embedded html file
            // Create the ResourceHandler. It is the object that will actually
            // handle the request for a given file. It is
            // a Jetty Handler object so it is suitable for chaining with other
            // handlers as you will see in other examples.
            ResourceHandler resource_handler = new ResourceHandler();
            // Configure the ResourceHandler. Setting the resource base
            // indicates where the files should be served out of.
            // In this example it is the current directory but it can be
            // configured to anything that the jvm has access to.
            resource_handler.setDirectoriesListed(true);
            resource_handler.setWelcomeFiles(new String[]{"index.html"});
            resource_handler.setResourceBase(".");

            // server.setHandler(new HandlerResource());
            ContextHandler context = new ContextHandler("/finder");
            // context.setContextPath("/");
            context.setHandler(new HandlerResource());
            ContextHandler contextMA = new ContextHandler("/mutilalign");
            contextMA.setHandler(new MutilAlignHandler());
            ContextHandler contextFS = new ContextHandler("/fakestruct");
            contextFS.setHandler(new FakeStructureHandler());

            // Add the ResourceHandler to the server.
            // ContextHandlerCollection contexts = new
            // ContextHandlerCollection();
            HandlerList contexts = new HandlerList();
            contexts.setHandlers(new Handler[]{resource_handler, context,
                    contextMA, contextFS});
            server.setHandler(contexts);
            server.start();
            server.join();

        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
}
