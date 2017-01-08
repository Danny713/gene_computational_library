/**
* GenEx GUI GraphFramer
* Frame containing SingleGraph/DoubleGraph and GraphOptionPanel
*
* Depends: SingleGraph/DoubleGraph, GraphOptionPanel, Java Swing
*
* This class was taken from the GenEx Project - written by Rachel Xu '17
**/

/*import libraries*/
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import edu.uci.ics.screencap.Dump;
import edu.uci.ics.screencap.PNGDump;


public class GraphFramer extends JFrame {
 /**
 * instance variables
 **/
 private JPanel glassPane;

 /**
 * constructor
 **/
 public GraphFramer(SingleGraph sg, double pval, double tau) {
  //create option panel
  //GraphOptionPanel gop = new GraphOptionPanel(sg, pval, tau, gd, this);

  //set frame attributes
  this.setLayout(new BorderLayout());
  this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

  //add components
  //this.add(gop, BorderLayout.LINE_START);
  this.add(sg, BorderLayout.CENTER);
  this.addGlassPane();


  this.pack();
  this.setVisible(true);
 }

 /**
 * add glass pane
 * copied from GenEx
 **/
 private void addGlassPane() {
  glassPane = new JPanel() {
   @Override
   public void paintComponent(Graphics g) {
    g.setColor(new Color(220, 220, 220, 200));
    g.fillRect(0, 0, this.getWidth(), this.getHeight());
   }
  };
  glassPane.setOpaque(false);
  glassPane.setLayout(new GridBagLayout());
  //wait label
  JLabel waitLabel = new JLabel("<html><div align = \"center\">Please wait.<br />Loading...</div></html>");
  waitLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 30));
  glassPane.add(waitLabel);

  //set glass pane
  this.setGlassPane(glassPane);
 }

 /**
 * glass pane visibility
 * copied from GenEx
 **/
 /*public void setVisibleGlass(boolean vis) {
  glassPane.setVisible(vis);
  this.repaint();
  this.revalidate();
 }*/
}