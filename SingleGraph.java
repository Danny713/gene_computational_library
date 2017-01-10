/**
 * This class was taken directly from the GenEx Project
 *
 * GenEx GUI SingleGraph
 * Creates a panel with the interactive graph, with the JUNG package
 * 
 * Depends: JUNG: UndirectedSparseGraph, VisualizationViewer + decorators, KKLayout, EditingModalGraphMouse, ResizeListener, Apache: Factory, Transformer
 **/

/*graph libraries*/
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.algorithms.layout.*;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.visualization.*;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;
import edu.uci.ics.jung.visualization.decorators.EdgeShape;
import edu.uci.ics.jung.visualization.control.*;
import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.Transformer;
/*GUI libraries*/
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.geom.*;
import java.awt.geom.Point2D;
import javax.swing.*;
import javax.swing.border.*;
import javax.imageio.ImageIO;
import java.io.*;
import java.awt.event.*;
/*other*/
import java.util.*;
import java.util.Collection;
import java.util.Vector;
import java.util.HashSet;
import edu.princeton.cs.introcs.*;


public class SingleGraph extends JPanel {
    //everything is protected to allow for extensibility (DoubleGraph)
    /**
     * constants
     **/
    //default dimensions
    public static final int DEFAULT_WIDTH = 700;
    public static final int DEFAULT_HEIGHT = 700;
    public static final int TOO_MANY = 100; //begin hiding labels
    public static final int REALLY_TOO_MANY = 1000; //disable vertex sizing and coloring scale
    protected final int MENUHEIGHT = 20;
    //colors
    protected final Color EDGE_COLOR = Color.GRAY;
    protected final int MAX_VERT = 20; //vertex size
    protected final int MIN_VERT = 5;
    protected final int MAX_RED = 255; //vertex color
    protected final int MIN_RED = 0;
    //other visual aspects
    protected final float EDGEWEIGHT = 2.5f; //maximum width of edge
    //protected final double EDGEWEIGHT = 2.5;
    
    /**
     * instance variables
     **/
    //graph construction
    protected KKLayout<String, Integer> l;
    protected Graph<String, Integer> g;
    protected Vector<String> isolates;
    protected Vector<Float> edgeWeights;
    //graph prettiness
    protected VisualizationViewer<String, Integer> vv;
    protected Transformer<Integer, Paint> edgePaint;
    protected Transformer<String, Paint> vertexPaint;
    protected Transformer<String, Shape> vertexShape;
    protected Transformer<String, String> vertexLabel;
    protected Transformer<Integer, Stroke> edgeStroke;
    //interactive functionality
    //protected EditingModalGraphMouse<String, Integer> gm;
    protected ResizeListener rl;
    //other variables
    protected String TITLE;
    protected int WIDTH;
    protected int HEIGHT;
    protected boolean needsKey; //whether or not to include a key
    protected boolean canToggle; //whether or not to include a toggle button
    protected int V; //largest possible number of vertices in graph
    protected HashSet<String> highlight;
    
    
    /**
     * constructor
     **/ 
    public SingleGraph() {
        ;
    }
    public SingleGraph(Correlate.EdgeList data, String title) {
        //check
        if (!data.isValid()) {
            throw new NullPointerException();
        }
        
        //initialize
        this.TITLE = title;
        highlight = new HashSet<String>();
        WIDTH = DEFAULT_WIDTH;
        HEIGHT = DEFAULT_HEIGHT;
        needsKey = false;
        canToggle = false;
        this.setLayout(new BorderLayout());
        
        createGraph(data); //create graph
        createLayout(); //create graph layout
        setAppearance(); //set graph appearance
        addMenuBar(); //add menu bar
        
        this.add(vv, BorderLayout.CENTER); //add viewer
        l.lock(false); //unlock vertices
        
        addResizeListener(); //add resizing listener
    }
    
    
    /**
     * create graph
     **/ 
    protected void createGraph(Correlate.EdgeList data) {
        g = new UndirectedSparseGraph<String, Integer>();
        V = data.names().length;
        
        addEdges(data); //add all edges & connected vertices
        fillIsolates(data); //get isolates
    }
    
    /*add edges to graph*/
    protected void addEdges(Correlate.EdgeList data) {
        edgeWeights = new Vector<Float>(V);
        
        for (int i = 0; i < data.E(); i++) {
            //get vertex names
            String first = data.getFirstName(i);
            String second = data.getSecondName(i);
            
            //add vertices -- does nothing if graph already contains vertices
            g.addVertex(first);
            g.addVertex(second);
            
            //add edge
            g.addEdge(i, first, second);
            edgeWeights.add(data.values()[i]);
        }
    }
    
    /*fill in isolates vector -
     case sensitive, but shouldn't be an issue*/
    protected void fillIsolates(Correlate.EdgeList data) {
        String[] names = data.names();
        
        //whoo fast lookup
        HashSet<String> verts = new HashSet<String>(g.getVertices());
        isolates = new Vector<String>();
        
        //add to isolates vector
        for (int i = 0; i < names.length; i++) {
            if (!verts.contains(names[i])) {
                isolates.add(names[i]);
            }
        }
    }
    
    
    /**
     * create graph layout
     **/
    public void createLayout() {
        Transformer<Integer, Double> t = new Transformer<Integer, Double>() {
            public Double transform(Integer edge) {
                return (1 - Math.pow(edgeWeights.elementAt(edge), 6));
            }
        };
        
        l = new KKLayout<String, Integer>(g, new DijkstraShortestPath<String, Integer>(g, t));
        //l = new KKLayout<String, Integer>(g);
        if (g.getVertexCount() > TOO_MANY) {
            l.setAdjustForGravity(false);
        }
        l.setSize(new Dimension(WIDTH, HEIGHT));
    }
    
    
    /**
     * set graph appearance
     **/
    protected void setAppearance() {
        //create edge/vertex/label/edge properties
        colorEdges();
        colorVertices();
        shapeVertices();
        edgeWidth();
        setLabels();
        
        //create & set visualization viewer
        vv = new VisualizationViewer<String, Integer>(l);
        vv.setPreferredSize(new Dimension(WIDTH,HEIGHT)); //set size
        setVVProperties();
        
        //add mouse
        addMouse();
    }
    
    /*set visualizer propertes*/
    protected void setVVProperties() {
        vv.setOpaque(true);
        vv.setBackground(Color.WHITE);
        vv.getRenderContext().setVertexLabelTransformer(vertexLabel);
        vv.getRenderContext().setVertexShapeTransformer(vertexShape);
        vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint);
        vv.getRenderContext().setEdgeDrawPaintTransformer(edgePaint);
        vv.getRenderContext().setEdgeShapeTransformer(new EdgeShape.Line<String, Integer>());
        vv.getRenderContext().setEdgeStrokeTransformer(edgeStroke);
    }
    
    /*add mouse for picking/choosing/translating graph*/
    protected void addMouse() {
        //mouse vertex factory
        Factory<String> vertexFactory = new Factory<String>() {
            public String create() {
                String newString = "";
                return newString;
            }
        };
        
        //mouse edge factory
        Factory<Integer> edgeFactory = new Factory<Integer>() {
            public Integer create() {
                int x = 0;
                return x;
            }
        };
        
        //mouse creation
        //gm = new EditingModalGraphMouse<String, Integer>(vv.getRenderContext(), vertexFactory, edgeFactory);
        //vv.setGraphMouse(gm);
        //gm.setMode(ModalGraphMouse.Mode.TRANSFORMING);
    }
    
    
    /**
     * create graph visualizer properties
     **/
    /*color edges*/
    protected void colorEdges() {
        edgePaint = 
            new Transformer<Integer, Paint>() {
            public Paint transform(Integer edge) {
                Color color = EDGE_COLOR;
                return color;
            }
        };
    }
    
    /*color vertices*/
    protected void colorVertices() {
        //if too many vertices, default to black
        if (g.getVertexCount() > REALLY_TOO_MANY) {
            vertexPaint =
                new Transformer<String, Paint>() {
                public Paint transform(String vertex) {
                    return Color.BLACK;
                }
            };
        }
        //otherwise
        else {
            vertexPaint =
                new Transformer<String, Paint>() {
                public Paint transform(String vertex) {
                    int number = g.getIncidentEdges(vertex).size() ;
                    int red = (int) MIN_RED;
                    int verts = g.getVertexCount();
                    if (g.getVertexCount() <= 0) {
                        verts = 1;
                    }
                    float step = (float) ((MAX_RED - MIN_RED) / (float) verts);
                    float scale = step * number + MIN_RED;
                    red = (int) (MIN_RED + step * number);
                    if (red > MAX_RED) red = (int) MAX_RED;
                    Color color = new Color((int) red, 0, 0);
                    
                    return color;
                }
            };
        }
    }
    
    /*shape vertices*/
    protected void shapeVertices() {
        //if too many vertices, default to minimum size
        if (g.getVertexCount() > REALLY_TOO_MANY) {
            vertexShape = 
                new Transformer<String, Shape>() {
                public Shape transform(String vertex) {
                    Ellipse2D.Float shape = new Ellipse2D.Float(-MIN_VERT/2, -MIN_VERT/2, MIN_VERT, MIN_VERT);
                    return shape;
                }
            };
        }
        //otherwise, go for scaling
        else {
            vertexShape = 
                new Transformer<String, Shape>() {
                public Shape transform(String vertex) {
                    
                    int number = g.getIncidentEdges(vertex).size();
                    int verts = g.getVertexCount();
                    if (g.getVertexCount() <= 0) {
                        verts = 1;
                    }
                    float step = ((float) (MAX_VERT - MIN_VERT) / (float) verts);
                    float size = step * number + MIN_VERT;
                    if (size > MAX_VERT) {
                        size = MAX_VERT;
                    }
                    Ellipse2D.Float shape = new Ellipse2D.Float(-size/2, -size/2, size, size);
                    return shape;
                }
            };
        }
    }
    
    /*add labels*/
    protected void setLabels() {
        if (g.getVertexCount() < TOO_MANY) {
            vertexLabel = new ToStringLabeller<String>();
        }
        else {
            vertexLabel = 
                new Transformer<String, String>() {
                public String transform(String vertex) {
                    return "";
                }
            };
        }
    }
    
    /*hide labels*/
    public void hideLabels() {
        vertexLabel = 
            new Transformer<String, String>() {
            public String transform(String vertex) {
                return "";
            }
        };
        vv.getRenderContext().setVertexLabelTransformer(vertexLabel);
        this.repaint();
    }
    
    /*show labels*/
    public void showLabels() {
        vertexLabel = new ToStringLabeller<String>();
        vv.getRenderContext().setVertexLabelTransformer(vertexLabel);
        this.repaint();
    }
    
    /*set edge stroke width*/
    protected void edgeWidth() {
        edgeStroke = 
            new Transformer<Integer, Stroke>() {
            public Stroke transform(Integer e) { 
                Float w = edgeWeights.elementAt(e);
                return new BasicStroke((EDGEWEIGHT * w * w * w * w));
            }
        };
    }
    
    
    /**
     * add menu bar
     **/
    protected void addMenuBar() {
        /*JMenuBar menuBar = new JMenuBar();
        JMenu modeMenu = gm.getModeMenu();
        modeMenu.setText("Select Mode of Graph Interaction");
        modeMenu.setIcon(null);
        modeMenu.setPreferredSize(new Dimension(WIDTH,MENUHEIGHT)); 
        menuBar.add(modeMenu);
        this.add(menuBar, BorderLayout.PAGE_START);*/
    }
    
    
    /**
     * add resize listener
     **/
    //add resize listener
    protected void addResizeListener() {
        rl = new ResizeListener(this);
        this.addComponentListener(rl);
    }
    
    
    /**
     * build a key panel (but not actually!)
     **/
    public JPanel keyPanel() {
        System.out.println("SingleGraph does not needs a key panel!");
        return null;
    }
    /**
     * toggle colors? (but not actually!)
     **/
    public void toggleColors() {
        ;
    }
    
    
    /**
     * access methods
     **/
    public KKLayout<String, Integer> getGraphLayout() {
        return l;
    }
    public VisualizationViewer<String, Integer> getViewer() {
        return vv;
    }
    public Vector<String> getIsolates() {
        return isolates;
    }
    public int vertexCount() {
        return g.getVertexCount();
    }
    public int edgeCount() {
        return g.getEdgeCount();
    }
    public String getTitle() {
        return TITLE;
    }
    public boolean needsKey() {
        return needsKey;
    }
    public boolean canToggle() {
        return canToggle;
    }
    public Collection<String> getVertices() {
        return g.getVertices();
    }
    
    
    /**
     * highlight certain vertices
     **/
    /*public void highlightVertices(DefaultListModel<String> toHighlight) {
        highlight.clear();
        for (Enumeration<String> e = toHighlight.elements(); e.hasMoreElements();) {
            String s = e.nextElement();
            highlight.add(s);
        }
        vertexPaint =
            new Transformer<String, Paint>() {
            public Paint transform(String vertex) {
                Color color = Color.BLACK;
                if (highlight.contains(vertex)) {
                    color = Color.GREEN;
                }
                return color;
            }
        };
        vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint);
        this.repaint();
    }
    /*stop highlight vertices*/
    public void stopHighlighting() {
        highlight.clear();
        colorVertices();
        vv.getRenderContext().setVertexFillPaintTransformer(vertexPaint);
        this.repaint();
    }
}



/**
 * resize listener
 **/
class ResizeListener extends ComponentAdapter {
    SingleGraph sg;
    
    public ResizeListener(SingleGraph sg) {
        this.sg = sg;
    }
    
    public void componentResized(ComponentEvent e) {
        if ((sg.getGraphLayout() != null) && (sg.getViewer() != null)) {
            Dimension d = new Dimension(sg.getViewer().getWidth(), sg.getViewer().getHeight());
            sg.getGraphLayout().setSize(d);
            sg.getViewer().setSize(d);
            sg.repaint();
        }
    }
}