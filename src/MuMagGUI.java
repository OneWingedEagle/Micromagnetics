
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.log10;
import static java.lang.System.out;
import javax.swing.*;


import java.awt.event.*;
import java.awt.*;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Scanner;

public class MuMagGUI extends JFrame{

	public		JPanel		ParaTab,DiscretizationTab,ModificationsTab,FieldandInitialmTab,RunTab,pp1;
	public JTextArea textArea=new JTextArea();
	public myJLabel lbK1,lbK2;
	public  myJTextField tfMs,tfA,tfK1,tfK2,tfgamma,tfalpha,tfdt,tfbetta;
	public myJTextField 	tfThetaMean,tfThetaRange,tfPhiMean,tfPhiRange,tfThetaCx,tfThetaCd,tfIterMax,tfErrorMax;
	public  myJTextField tfLx,tfLy,tfLz,tfNx,tfNy,tfNz,tdx,tdy,tdz;
	public  myJTextField tfRootDirectory,tfMeshFile,tfmagFile,tfPrevmagFile,tfHeffFile,tfHeffSurfFile,tfmagSurfFile,tfParametersFile;
	
	public	 myJTextField tfHa,tfThetaHa,tfPhiHa;
	public int Nx=1,Ny=1,Nz=1;
	public double Ms,A,K1,K2,gamma,alpha,betta,dt,Lx,Ly,Lz,dx,dy,dz;
	public double Ha,PhiHa,ThetaHa,ThetaMean,ThetaRange,PhiMean,PhiRange,ThetaCx,ThetaCd;
	public  myJButton Discretize,Run;
	public JPanel pp3,pp4;
	public double lx,ErrorMax;
	public double mu0=Math.PI*4e-7;
	public boolean disc;
	JComboBox dynSelector,aniSelector;
	public boolean Rot=true,LLG,PreviousInitial,ModifEex,ModifEd,CubicAni=true,prevPara;
	public JCheckBox chkPrev,chkModifEex,chkModifEd,chkPrevPara;
	public int IterMax;
	Vect uniAniAxis=new Vect(3);
		
	 
	 MuMagGUI() {
		 
		 JTabbedPane tbPanel = new JTabbedPane();
			getContentPane().add(tbPanel);
			ParaTab = new JPanel(new GridLayout(3,1,10,10));
			DiscretizationTab = new JPanel(new GridLayout(3,1,10,10));
			FieldandInitialmTab = new JPanel(new GridLayout(3,1,10,10));
			ModificationsTab = new JPanel(new GridLayout(5,1,10,10));
			RunTab = new JPanel(new GridLayout(1,2,10,10));
			
			
	 
			tbPanel.setFont(new Font("Arial", 1, 12));
			tbPanel.addTab("Parameters",ParaTab);
			tbPanel.addTab("Dircretization", DiscretizationTab);
			tbPanel.addTab("Modifications", ModificationsTab);
			tbPanel.addTab("Applied Field and Initial values of m", FieldandInitialmTab);
			tbPanel.addTab("Run and Write", RunTab);
			
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			setTitle("Micromagnetic Simulation"); 
			setSize(850,650);
			setLocation(10, 10);

		 
		    
	 //================================================================ redirecting console to text area
		
			textArea.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			textArea.setEditable(false);;
			textArea.setBorder(BorderFactory.createLineBorder(Color.blue,1));
			   JScrollPane scrollPane = new JScrollPane(textArea);
			   scrollPane.setBorder(BorderFactory.createCompoundBorder(
						BorderFactory.createTitledBorder("Progress"),
						BorderFactory.createEmptyBorder(10,5,5,5)));
			   scrollPane.setPreferredSize(new Dimension(400,450));

			
	//================================================================ Creating ParaTab
			
			tfMs=new myJTextField();
			tfA=new myJTextField();
			tfK1=new myJTextField();
			tfK2=new myJTextField();
			tfgamma=new myJTextField();
			tfalpha=new myJTextField();
			tfdt=new myJTextField();
			tfbetta=new myJTextField();
			tfLx=new myJTextField();
			tfLy=new myJTextField();
			tfLz=new myJTextField();

			myJLabel lbMs=new  myJLabel(" Saturation Magnetization: Ms"   , myJLabel.RIGHT);
			myJLabel lbMsunit=new  myJLabel("A/m"   , myJLabel.LEFT);

			myJLabel lbA=new myJLabel("Exchange Stiffness Constant: A", myJLabel.RIGHT);
			myJLabel lbAunit=new  myJLabel("J/m"   , myJLabel.LEFT);

			lbK1=new myJLabel("First Cubic Anisotropy Constant: K1", myJLabel.RIGHT);
			myJLabel lbK1unit=new  myJLabel("J/m3"   , myJLabel.LEFT);

			lbK2=new myJLabel("Second Cubic Anisotropy Constant: K2", myJLabel.RIGHT);
			myJLabel lbK2unit=new  myJLabel("J/m3"   , myJLabel.LEFT);

			myJTextField EasyAxes=new myJTextField();
			myJLabel lbgamma=new  myJLabel(" Gyromagnetic ratio: gamma"   , myJLabel.RIGHT);
			tfalpha=new myJTextField();
			tfgamma=new myJTextField();
			tfdt=new myJTextField();
			myJLabel lbgammaunit=new  myJLabel("m/As"   , myJLabel.LEFT);
			myJLabel lbalpha=new myJLabel("Damping Constant: alpha", myJLabel.RIGHT);
			myJLabel lbalphaunit=new  myJLabel(""   , myJLabel.LEFT);
			myJLabel lbdt=new myJLabel("TimeStep: dt", myJLabel.RIGHT);
			myJLabel lbdtunit=new  myJLabel("Picoseconds"   , myJLabel.LEFT);

			myJLabel lbbetta=new myJLabel("Relaxation Factor: 0<betta<1", myJLabel.RIGHT);

			pp1 = new JPanel(new GridLayout(6,3,10,10));
			JPanel pp2 = new JPanel(new GridLayout(4,3,10,10));
			pp3 = new JPanel(new GridLayout(4,3,10,10));
			pp4 = new JPanel(new GridLayout(4,3,10,10));

			ParaTab.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);

			tfMs.setText("1.72e6");
			tfA.setText("2.01e-11");
			tfK1.setText("4.8e4");
			tfK2.setText("1.5e3");
			tfgamma.setText("2.21e5");
			tfalpha.setText("0.2");
			tfdt.setText("0.2");
			tfbetta.setText("0.2");
			chkPrevPara=new JCheckBox("Load Previous Parameters");
			chkPrevPara.addItemListener(new  ItemListener(){
				public void itemStateChanged(ItemEvent e){
					if(chkPrevPara.isSelected()){
						prevPara=true;
						loadParameters();
											}
					else
						prevPara=false;
					

				}

			});
			
			myJLabel lbAnisotropy=new myJLabel("Anisotropy type", myJLabel.RIGHT);
			String[] anisotropy = { " Cubic( x,y, and z axes)"," Uniaxial( x axis)"," Uniaxial( y axis)"," Uniaxial( z axis)"};
			 aniSelector =new JComboBox(anisotropy);
			aniSelector.setBackground(new Color(215, 255, 225));
			aniSelector.setSelectedIndex(0);
			aniSelector.addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent e)
				{				
					int Selection=aniSelector.getSelectedIndex();
					if(Selection==1){
						CubicAni=false;
						uniAniAxis=new Vect(3);
						uniAniAxis.el[0]=1;
					lbK1.setText(" Firs Uniaxial Anisotropy Constrant");
					lbK2.setText(" Second Uniaxial Anisotropy Constrant");
								}

					else if(Selection==2){
						CubicAni=false;
						uniAniAxis=new Vect(3);
						uniAniAxis.el[1]=1;
						lbK1.setText(" Firs Uniaxial Anisotropy Constrant");
						lbK2.setText(" Second Uniaxial Anisotropy Constrant");
		
					}
					else if(Selection==3){
						CubicAni=false;
						uniAniAxis=new Vect(3);
						uniAniAxis.el[2]=1;
						lbK1.setText(" Firs Uniaxial Anisotropy Constrant");
						lbK2.setText(" Second Uniaxial Anisotropy Constrant");
		
					}
					else {
						CubicAni=true;
						lbK1.setText(" Firs Cubic Anisotropy Constrant");
						lbK2.setText(" Second Cubic Anisotropy Constrant");
					}
				}
			});
			EasyAxes.setText("x, y, and z axes");
			EasyAxes.setEditable(false);
			myJLabel lblexch=new myJLabel("Exchange length", myJLabel.RIGHT);
			final myJTextField lexch=new myJTextField();
			myJLabel lblexUnit=new myJLabel("nm", myJLabel.LEFT);
			myJButton CalExchlength=new myJButton("Calculate exchange length");
			
			CalExchlength.addActionListener(new  ActionListener(){
				public void actionPerformed(ActionEvent e){
					Ms=Double.parseDouble(tfMs.getText());  
					A=Double.parseDouble(tfA.getText()); 
					lx=Math.sqrt(2*A/(mu0*Ms*Ms));			
					DecimalFormat f = new DecimalFormat("0.00");
					double lxn=lx*1e9;
					lexch.setText(f.format(lxn));
				}
			});

			myJLabel lbDynamic=new myJLabel("Dynamic: Rotation or LLG?", myJLabel.RIGHT);
			String[] dynamic = { "  Rotation  ", "  LLG"};

			dynSelector =new JComboBox(dynamic);
			dynSelector.setBackground(new Color(215, 255, 225));
			dynSelector.setFont(new Font("Arial", 1, 15));
			dynSelector.setSelectedIndex(0);
			dynSelector.addItemListener(new ItemListener() {
				public void itemStateChanged(ItemEvent e)
				{				
					int Selection=dynSelector.getSelectedIndex();
					if(Selection==1){
						LLG=true;
						Rot=false;
						ParaTab.remove(pp3);;
						ParaTab.add(pp4);;
						ParaTab.revalidate();
						ParaTab.repaint();
					}

					else if(Selection==0){
						LLG=false;
						Rot=true;
						ParaTab.remove(pp4);;
						ParaTab.add(pp3);;
						ParaTab.revalidate();
						ParaTab.repaint();
					}
				}
			});
			pp1.add(lbMs);
			pp1.add(tfMs);
			pp1.add(lbMsunit);
			pp1.add(lbA);
			pp1.add(tfA);
			pp1.add(lbAunit);
			pp1.add(lbAnisotropy);
			pp1.add(aniSelector);
			pp1.add(new myJLabel());
			pp1.add(lbK1);
			pp1.add(tfK1);
			pp1.add(lbK1unit);
			pp1.add(lbK2);
			pp1.add(tfK2);
			pp1.add(lbK2unit);
			
			pp2.add(new myJLabel());
			pp2.add(chkPrevPara);
			pp2.add(CalExchlength);
			pp2.add(lblexch);
			pp2.add(lexch);
			pp2.add(lblexUnit);
			pp2.add(new myJLabel());
			pp2.add(new myJLabel());
			pp2.add(new myJLabel());
			pp2.add(lbDynamic);
			pp2.add(dynSelector);
			pp2.add(new myJLabel());


			pp3.add(lbbetta);
			pp3.add(tfbetta);
			pp3.add(new myJLabel());
			pp3.add(new myJLabel());
			pp3.add(new myJLabel());
			pp3.add(new myJLabel());
			pp3.add(new myJLabel());
			pp3.add(new myJLabel());
			pp3.add(new myJLabel());
			
			pp4.add(lbgamma);
			pp4.add(tfgamma);
			pp4.add(lbgammaunit);
			pp4.add(lbalpha);
			pp4.add(tfalpha);
			pp4.add(lbalphaunit);
			pp4.add(lbdt);
			pp4.add(tfdt);
			pp4.add(lbdtunit);
			pp4.add(new myJLabel());
			pp4.add(new myJLabel());
			pp4.add(new myJLabel());


			ParaTab.add(pp1);
			ParaTab.add(pp2);
			ParaTab.add(pp3);

		//================================================================ Creating DiscretizationTab
			
			JPanel gp1 = new JPanel(new GridLayout(4,4,10,10));
			JPanel gp2 = new JPanel(new GridLayout(4,4,10,10));
			DiscretizationTab.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
			DiscretizationTab.add(gp1);
			DiscretizationTab.add(gp2);
			DiscretizationTab.add(new JPanel());

			tfLx.setText("40.0");
			tfLy.setText("40.0");
			tfLz.setText("5.0");
			
			myJLabel lbLx=new  myJLabel("Length in x direction (nm)"   , myJLabel.RIGHT);
			myJLabel lbLy=new  myJLabel("Length in y direction (nm)"   , myJLabel.RIGHT);
			myJLabel lbLz=new  myJLabel("Length in z direction (nm)"   , myJLabel.RIGHT);
			myJLabel lbdx=new  myJLabel("Cell size along x: dx (nm)"   , myJLabel.RIGHT);
			myJLabel lbdy=new  myJLabel("Cell size along y: dy (nm)"   , myJLabel.RIGHT);
			myJLabel lbdz=new  myJLabel("Cell size along z: dz (nm)"   , myJLabel.RIGHT);
			myJLabel lbNx=new  myJLabel("Divisons in x direction: Nx"   , myJLabel.RIGHT);
			myJLabel lbNy=new  myJLabel("Divisons in y direction: Ny"   , myJLabel.RIGHT);
			myJLabel lbNz=new  myJLabel("Divisons in z direction: Nz"   , myJLabel.RIGHT);

			tdx=new myJTextField();
			tdy=new myJTextField();
			tdz=new myJTextField();
		

			tfNx=new myJTextField();
			tfNy=new myJTextField();
			tfNz=new myJTextField();
			tfNx.setText("8");
			tfNy.setText("8");
			tfNz.setText("1");

			Discretize=new myJButton("Discretize!");
			
			Discretize.addActionListener(new  ActionListener(){
				public void actionPerformed(ActionEvent e){
				
					
					Lx=Double.parseDouble(tfLx.getText());  
					Ly=Double.parseDouble(tfLy.getText()); 
					Lz=Double.parseDouble(tfLz.getText());
					Nx=Integer.parseInt(tfNx.getText()); 
					Ny=Integer.parseInt(tfNy.getText()); 
					Nz=Integer.parseInt(tfNz.getText()); 
								
					if(Nx<1 || Ny<1 || Nz<1) {
						String msg="divisions can not be less than 1";
						JOptionPane.showMessageDialog(null, msg," ", JOptionPane. WARNING_MESSAGE);
					}
					else
					{
						dx=Lx/Nx;
						dy=Ly/Ny;
						dz=Lz/Nz;	
						
					DecimalFormat f = new DecimalFormat("0.00");
				
					tdx.setText(f.format(dx));
					tdy.setText(f.format(dy));
					tdz.setText(f.format(dz));	
					
					disc=true;
				}
				}
			});

			gp1.add(lbLx);
			gp1.add(tfLx);
			gp1.add(lbNx);
			gp1.add(tfNx);
			gp1.add(lbLy);
			gp1.add(tfLy);
			gp1.add(lbNy);
			gp1.add(tfNy);
			gp1.add(lbLz);
			gp1.add(tfLz);
			gp1.add(lbNz);
			gp1.add(tfNz);
			gp1.add(new myJLabel());
			gp1.add(Discretize);
			gp1.add(new myJLabel());
			gp1.add(new myJLabel());
			
			gp2.add(new myJLabel());
			gp2.add(lbdx);
			gp2.add(tdx);
			gp2.add(new myJLabel());
			gp2.add(new myJLabel());
			gp2.add(lbdy);
			gp2.add(tdy);
			gp2.add(new myJLabel());
			gp2.add(new myJLabel());
			gp2.add(lbdz);
			gp2.add(tdz);
			gp2.add(new myJLabel());
			gp2.add(new myJLabel());
			gp2.add(new myJLabel());
			gp2.add(new myJLabel());
		
	//================================================================ Creating ModificationsTab
			JPanel mp1 = new JPanel(new GridLayout(1,2,10,10));
			JPanel mp2 = new JPanel(new GridLayout(2,3,10,10));
			
			chkModifEex=new JCheckBox("Use modifed representation of excange energy",false);
			chkModifEd=new JCheckBox("Use modifed representation of demagnetizing energy",false);
			myJLabel lbThetaCx=new myJLabel("Critical angle for Exch. representation", myJLabel.RIGHT);
			myJLabel lbThetaCd=new myJLabel("Critical angle for Demag. representation", myJLabel.RIGHT);
			myJLabel lbThetaCxUnit=new myJLabel("degree", myJLabel.LEFT);
			myJLabel lbThetaCdUnit=new myJLabel("degree", myJLabel.LEFT);
			tfThetaCx=new myJTextField("180.0");
			tfThetaCd=new myJTextField("180.0");
			tfThetaCx.setEditable(false);
			tfThetaCd.setEditable(false);
			
			chkModifEex.setFont(new Font("Arial", 1, 15));
			chkModifEd.setFont(new Font("Arial", 1, 15));
			
			chkModifEex.addItemListener(new  ItemListener(){
				public void itemStateChanged(ItemEvent e){
					if(chkModifEex.isSelected()){
						ModifEex=true;
						tfThetaCx.setEditable(true);
						}
					else{
						tfThetaCx.setEditable(false);
						ModifEex=false;
					}

				}

			});
				
			chkModifEd.addItemListener(new  ItemListener(){
				public void itemStateChanged(ItemEvent e){
					if(chkModifEd.isSelected()){
						ModifEd=true;
						tfThetaCd.setEditable(true);
						}
					else{
						ModifEd=false;
						tfThetaCd.setEditable(false);
					}

				}

			});
			mp1.add(chkModifEex);
			mp1.add(chkModifEd);
			mp2.add(lbThetaCx);
			mp2.add(tfThetaCx);
			mp2.add(lbThetaCxUnit);
			mp2.add(lbThetaCd);
			mp2.add(tfThetaCd);
			mp2.add(lbThetaCdUnit);
			
			ModificationsTab.add(mp1);
			ModificationsTab.add(mp2);


	//================================================================ Creating FieldandInitialmTab
			
			JPanel fp1 = new JPanel(new GridLayout(6,2,10,10));
			JPanel fp2 = new JPanel(new GridLayout(1,1,10,10));
			JPanel fp3 = new JPanel(new GridLayout(6,3,10,10));

			FieldandInitialmTab.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
			FieldandInitialmTab.add(fp1);
			FieldandInitialmTab.add(fp2);
			FieldandInitialmTab.add(fp3);

		
			myJLabel lbHa=new myJLabel("Magnitude of the applied field", myJLabel.RIGHT);
			myJLabel lbThetaHa=new myJLabel("Polar angle of the applied field ", myJLabel.RIGHT);
			myJLabel lbPhiHa=new myJLabel("Azimuthal angle of the applied field", myJLabel.RIGHT);
			myJLabel lbHaUnit=new myJLabel("A/m", myJLabel.LEFT);
			myJLabel lbThetaHaUnit=new myJLabel("degree", myJLabel.LEFT);
			myJLabel lbPhiHaUnit=new myJLabel("degree", myJLabel.LEFT);


			
			tfHa=new myJTextField("0.0");
			tfThetaHa=new myJTextField("0.0");
			tfPhiHa=new myJTextField("0.0");
			myJLabel lbApplied=new myJLabel("Applied Field:");
			lbApplied.setFont(new Font("Times New Roman",1,16));
/*			lbApplied.setOpaque(true);
			lbApplied.setBackground(Color.green);*/
			
			fp1.add(new JLabel());
			fp1.add(new JLabel());
			fp1.add(new JLabel());
			fp1.add(lbApplied);
			fp1.add(new JLabel());
			fp1.add(new JLabel());
			fp1.add(lbHa);
			fp1.add(tfHa);
			fp1.add(lbHaUnit);
			fp1.add(lbThetaHa);
			fp1.add(tfThetaHa);
			fp1.add(lbThetaHaUnit);
			fp1.add(lbPhiHa);
			fp1.add(tfPhiHa);
			fp1.add(lbPhiHaUnit);
			fp1.add(new JLabel());
			fp1.add(new JLabel());
			fp1.add(new JLabel());
			
			myJLabel lbInitial=new myJLabel("Setting initial values:");
			lbInitial.setFont(new Font("Times New Roman",1,16));
			//bInitial=new myJButton("Initial value of m");

			myJLabel lbThetaMean=new myJLabel("Theta_mrean", myJLabel.RIGHT);
			myJLabel lbThetaRange=new myJLabel("Theta_random_range", myJLabel.RIGHT);
			lbThetaRange.setFont(new Font("Arial", 0, 14));
			myJLabel lbPhiMean=new myJLabel("Phi_mean", myJLabel.RIGHT);
			myJLabel lbPhiRange=new myJLabel("Phi_random_range", myJLabel.RIGHT);
			lbPhiRange.setFont(new Font("Arial", 0, 14));
			myJLabel lbThetaMeanUnit=new myJLabel("degree", myJLabel.LEFT);
			myJLabel lbThetaRangeUnit=new myJLabel("degree ", myJLabel.LEFT);
			myJLabel lbPhiMeanUnit=new myJLabel("degree", myJLabel.LEFT);
			myJLabel lbPhiRangeUnit=new myJLabel("degree", myJLabel.LEFT);

			tfThetaMean=new myJTextField("90.0");
			tfThetaRange=new myJTextField("0.0");
			tfPhiMean=new myJTextField("90.0");
			tfPhiRange=new myJTextField("0.0");
			/*bInitial.setHorizontalAlignment(0);
			bInitial.setFont(new Font("Times New Roman", 0, 20));*/
			chkPrev=new JCheckBox("Use the result of the previous run as initial valuse of m",false);
		
			chkPrev.setFont(new Font("Arial", 1, 15));
			chkPrev.addItemListener(new  ItemListener(){
				public void itemStateChanged(ItemEvent e){
					if(chkPrev.isSelected()){
						PreviousInitial=true;
						//bInitial.setEnabled(false);

						tfThetaMean.setEditable(false);
						tfThetaRange.setEditable(false);
						tfPhiMean.setEditable(false);
						tfPhiRange.setEditable(false);
					}
					else{
						PreviousInitial=false;
						//bInitial.setEnabled(true);
						tfThetaMean.setEditable(true);
						tfThetaRange.setEditable(true);
						tfPhiMean.setEditable(true);
						tfPhiRange.setEditable(true);
					}

				}

			});

			fp2.add(chkPrev);
			
			fp3.add(lbInitial);
			fp3.add(new JLabel());
			fp3.add(new JLabel());
		
			fp3.add(lbThetaMean);
			fp3.add(tfThetaMean);
			fp3.add(lbThetaMeanUnit);
			fp3.add(lbThetaRange);
			fp3.add(tfThetaRange);
			fp3.add(lbThetaRangeUnit);
			
			fp3.add(lbPhiMean);
			fp3.add(tfPhiMean);
			fp3.add(lbPhiMeanUnit);
			fp3.add(lbPhiRange);
			fp3.add(tfPhiRange);
			fp3.add(lbPhiRangeUnit);
			fp3.add(new JLabel());
			fp3.add(new JLabel());
			fp3.add(new JLabel());

	//================================================================ Creating RunTab		
			 
				myJLabel lbsaveMesh=new  myJLabel("Save mesh to"   , myJLabel.RIGHT);
				myJLabel lbRootDirectory=new  myJLabel("Root directory"   , myJLabel.RIGHT);
				myJLabel lbsavemag=new  myJLabel("Save m to"   , myJLabel.RIGHT);
				myJLabel lbsavem4ini=new  myJLabel("Save m for the next run to"   , myJLabel.RIGHT);
				myJLabel lbsaveHeff=new  myJLabel("Save Heff to"   , myJLabel.RIGHT);
				myJLabel lbsaveHeffSurf=new  myJLabel("Save Heff of sufrace to"   , myJLabel.RIGHT);
				myJLabel lbsavemSurf=new  myJLabel("Save m of sufrace to"   , myJLabel.RIGHT);
				myJLabel lbErrorMax=new  myJLabel("Error criterion( degree)"   , myJLabel.RIGHT);
				myJLabel lbIterMax=new  myJLabel("Max. iteration"   , myJLabel.RIGHT);
				
				tfRootDirectory=new myJTextField(System.getProperty("user.dir"));
				tfMeshFile=new myJTextField("mesh.txt");
				tfmagFile=new myJTextField("mag.txt");
				tfPrevmagFile=new myJTextField("previousm.txt");
				tfPrevmagFile.setEditable(false);
				tfHeffFile=new myJTextField("Heff.txt");
				tfHeffSurfFile=new myJTextField("HeffSurf.txt");
				tfmagSurfFile=new myJTextField("magSurf.txt");
				tfParametersFile=new myJTextField("parameters.txt");
				tfParametersFile.setEditable(false);
				tfIterMax=new myJTextField("100");
				tfErrorMax=new myJTextField("0.1");
				JPanel rp1 = new JPanel(new GridLayout(14,2,10,10));
				JPanel rp21 = new JPanel(new GridLayout(1,2,10,10));
				JPanel rp22 = new JPanel(new FlowLayout(0,5,5));
				 rp22.setBorder(BorderFactory.createEmptyBorder(50,5,5,0));
				
				RunTab.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
				RunTab.add(rp21);
				myJLabel lbsaveTo=new myJLabel("Save output to :", myJLabel.RIGHT);
				lbsaveTo.setFont(new Font("Arial", 1, 16));
				
				 Run=new myJButton("Run");
				Run.setPreferredSize(new Dimension(100,40));
				 

				
				
				rp1.add(lbIterMax);
				rp1.add(tfIterMax);
				rp1.add(lbErrorMax);
				rp1.add(tfErrorMax);
				rp1.add(lbsaveTo);	
				rp1.add(new myJLabel());
				rp1.add(lbRootDirectory);
				rp1.add(tfRootDirectory);
				rp1.add(lbsaveMesh);
				rp1.add(tfMeshFile);
				rp1.add(lbsavemag);
				rp1.add(tfmagFile);
				rp1.add(lbsavem4ini);
				rp1.add(tfPrevmagFile);
				rp1.add(lbsaveHeff);
				rp1.add(tfHeffFile);
				rp1.add(lbsaveHeffSurf);
				rp1.add(tfHeffSurfFile);
				rp1.add(lbsavemSurf);
				rp1.add(tfmagSurfFile);
				rp1.add(new myJLabel("Save Parameters to",myJLabel.RIGHT));
				rp1.add(tfParametersFile);
				rp1.add(new myJLabel());
				rp1.add(new myJLabel());
				rp1.add(new myJLabel());
				rp1.add(new myJLabel());
				rp1.add(new myJLabel());
				rp1.add(new myJLabel());

				rp22.add(Run);
				rp22.add(new JLabel());
				rp22.add(scrollPane);
				
				rp21.add(rp1);
				rp21.add(rp22);

	}



	
	public void loadParameters(){
		String folder=System.getProperty("user.dir");
		try{
			FileReader fr = new FileReader(folder+"\\parameters.txt");
			Scanner scr = new Scanner(fr);
			while(!scr.next().equals("A:")){}
			tfA.setText(scr.next());
			while(!scr.next().equals("Ms:")){}
			tfMs.setText(scr.next());
			while(!scr.next().equals("Anisotropy:")){}
			if(scr.next().equals("(x,y,"))
				aniSelector.setSelectedIndex(0);
		
			else if(scr.next().equals("x"))
					aniSelector.setSelectedIndex(1);
			else if(scr.next().equals("y"))
				aniSelector.setSelectedIndex(2);
			else
				aniSelector.setSelectedIndex(3);
			while(!scr.next().equals("K1:")){}
			tfK1.setText(scr.next());
			while(!scr.next().equals("K2:")){}
			tfK2.setText(scr.next());
		
			while(!scr.next().equals("Dynamic:")){}
			if(scr.next().equals("LLG"))
				dynSelector.setSelectedIndex(1);
			else
				dynSelector.setSelectedIndex(0);

			
			scr.close();
			fr.close();
			}
			catch(IOException e){System.err.println("Error in reading previous parameters");
			}
	}
}




