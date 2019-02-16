import static java.lang.Math.PI;
import static java.lang.Math.log10;
import static java.lang.Math.pow;
import static java.lang.System.out;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Scanner;
import javax.swing.JOptionPane;
import javax.swing.JTextArea;

public class Main{
	private String  rootDirectory,meshFile,magFile,prevmagFile,HeffFile,HeffSurfFile,magSurfFile,parametersFile;
	private MuMagGUI gui;
	private Model model;
	private int iterMax,runCounter=0;
	private double Ha,thetaHa,phiHa; 
	private double mu0=4*PI*1e-7;
	private double errorMax;
	private double thetaMean,thetaRandomRange,phiMean,phiRandomRange,thetaCx,thetaCd;

	Main(){	  

		gui=new MuMagGUI();
		model=new Model();
		gui.setVisible(true);

		gui.Run.addActionListener(new  ActionListener(){
			public void actionPerformed(ActionEvent e){

				if(!gui.PreviousInitial && !gui.disc){
					String msg="Model has not been discretized yet. Please open the Discritization tab and descretize";
					JOptionPane.showMessageDialog(null, msg," ", JOptionPane. WARNING_MESSAGE);
				}
				else if(!gui.LLG && gui.betta>1){
					model.betta=Double.parseDouble(gui.tfbetta.getText());
					String msg="Relaxation factor must not be greater than one";
					JOptionPane.showMessageDialog(null, msg," ", JOptionPane. WARNING_MESSAGE);
				}
				else{
					runCounter++;

					model.Ms=Double.parseDouble(gui.tfMs.getText());
					model.mu0Ms=mu0*model.Ms;
					model.A=Double.parseDouble(gui.tfA.getText());
					model.K1=Double.parseDouble(gui.tfK1.getText());
					model.K2=Double.parseDouble(gui.tfK2.getText());
					model.gamma=Double.parseDouble(gui.tfgamma.getText());
					model.alpha=Double.parseDouble(gui.tfalpha.getText());
					model.dt=Double.parseDouble(gui.tfdt.getText())*1e-12;
					model.betta=Double.parseDouble(gui.tfbetta.getText());
					model.uniAniAxis.el=Arrays.copyOf(gui.uniAniAxis.el, 3);
					Ha=Double.parseDouble(gui.tfHa.getText());
					phiHa=Double.parseDouble(gui.tfPhiHa.getText());
					thetaHa=Double.parseDouble(gui.tfThetaHa.getText());	
					iterMax=Integer.parseInt(gui.tfIterMax.getText());
					errorMax=Double.parseDouble(gui.tfErrorMax.getText());
					rootDirectory=gui.tfRootDirectory.getText();
					meshFile=rootDirectory+"\\"+gui.tfMeshFile.getText();
					magFile=rootDirectory+"\\"+gui.tfmagFile.getText();
					prevmagFile=rootDirectory+"\\"+gui.tfPrevmagFile.getText();
					magSurfFile=rootDirectory+"\\"+gui.tfmagSurfFile.getText();
					HeffFile=rootDirectory+"\\"+gui.tfHeffFile.getText();
					HeffSurfFile=rootDirectory+"\\"+gui.tfHeffSurfFile.getText();
					parametersFile=rootDirectory+"\\"+gui.tfParametersFile.getText();

					if(!gui.PreviousInitial){
						thetaMean=Double.parseDouble(gui.tfThetaMean.getText());
						thetaRandomRange=Double.parseDouble(gui.tfThetaRange.getText());
						phiMean=Double.parseDouble(gui.tfPhiMean.getText());
						phiRandomRange=Double.parseDouble(gui.tfPhiRange.getText());
					}				
					if(gui.ModifEex) thetaCx=Double.parseDouble(gui.tfThetaCx.getText());
					if(gui.ModifEd) thetaCd=Double.parseDouble(gui.tfThetaCd.getText());
					model.CubicAni=gui.CubicAni;
					model.LLG=gui.LLG;
					model.PreviousInitial=gui.PreviousInitial;
					model.ModifEex=gui.ModifEex;
					model.ModifEd=gui.ModifEex;



					if(gui.PreviousInitial){
						try{
							FileReader fr = new FileReader(prevmagFile);
							Scanner scr = new Scanner(fr);
							model.Lx=scr.nextDouble();
							model.Ly=scr.nextDouble();
							model.Lz=scr.nextDouble();
							model.Ny=scr.nextInt();
							model.Nx=scr.nextInt();
							model.Nz=scr.nextInt();
							model. dxnm=model.Lx/model.Nx;
							model.dynm=model.Ly/model.Ny;
							model.dznm=model.Lz/model.Nz;

							scr.close();
							fr.close();
						}
						catch(IOException e2){System.err.println("Error in reading the previous result");
						}

					}
					else{
						model.Lx=Double.parseDouble(gui.tfLx.getText());
						model.Ly=Double.parseDouble(gui.tfLy.getText());
						model.Lz=Double.parseDouble(gui.tfLz.getText());
						model.Nx=Integer.parseInt(gui.tfNx.getText());
						model.Ny=Integer.parseInt(gui.tfNy.getText());
						model.Nz=Integer.parseInt(gui.tfNz.getText());
						model.dxnm=model.Lx/model.Nx;
						model.dynm=model.Ly/model.Ny;
						model.dznm=model.Lz/model.Nz;
					}
					model.I=model.Ny;
					model.J=model.Nx;
					model.K=model.Nz;

					Console.redirectOut(gui.textArea);
					writeParameters();

					Thread muMagThread=new Thread(){
						public void run(){
							model.setCells();
							double t1,t2,iterTime=0;
							out.println();
							int scale = 1;
							if(model.PreviousInitial)
								model.setInitialPrev(prevmagFile);
							else
								model.setInitialNew(thetaMean,thetaRandomRange,phiMean,phiRandomRange);
							model.setDemagTensorFFT();
							model.setHext(Ha,thetaHa,phiHa);

							int it;
							out.println();
							out.println("Energy minimization started.");
							out.println();
							out.format("%s\t%s\t%s\t%s\n", new Object[] {
									"  Iteration", "Etotal", "Error", "time per iteration (seconds)"
							});
							for(it=1;it<=iterMax;it++){
								if(it>1 && model.AngMax<errorMax/180*PI) break;
								t1=System.currentTimeMillis();
								model.calHd();
								model.calHexch();
								model.calHani();
								model.calHeff();
								model.calEext();
								model.calEtotal();
								model.movem();
								t2=System.currentTimeMillis();

								iterTime=(t2-t1)/1000;
								if(it==1 || it%10==0){
									if(model.Etotal==0)
										scale=0;
									else
										scale = -(int)log10(model.Etotal) + 1;
									System.out.format("% d\t%.3fE%d\t%.3f\t%.3f\n",it, pow(10, scale) * model.Etotal, -scale, (model.AngMax / PI) * 180,iterTime);
								}
							}

							if(it%10!=0) 
								System.out.format("%s%d\t%.3fE%d\t%.3f\t%.3f\n","  ",it, pow(10, scale) * model.Etotal, -scale, (model.AngMax / PI) * 180,iterTime);
							System.out.println();
							if(it==iterMax+1)
								System.out.println("  Iteration completed. ");
							else
								System.out.println("  Desired convergence achieved.");
							System.out.println();
							model.writeMesh(meshFile);
							System.out.println((new StringBuilder(" Writing ")).append(meshFile).append(" completed.").toString());
							System.out.println();
							model.writemag(magFile);
							System.out.println((new StringBuilder(" Writing ")).append(magFile).append(" completed.").toString());
							System.out.println();
							model.writemForIni(prevmagFile);
							System.out.println((new StringBuilder(" Writing ")).append(prevmagFile).append(" completed.").toString());
							System.out.println();
							model.writemagSurf(magSurfFile);
							System.out.println((new StringBuilder(" Writing ")).append(magSurfFile).append(" completed.").toString());
							System.out.println();
							model.writeHeff(HeffFile);
							System.out.println((new StringBuilder(" Writing ")).append(HeffFile).append(" completed.").toString());
							System.out.println();
							model.writeHeffSurf(HeffSurfFile);
							System.out.println((new StringBuilder(" Writing ")).append(HeffSurfFile).append(" completed.").toString());
							System.out.println();
							System.out.println(" Writing files completed. ");
							out.println("============================================");
							out.println("============================================");

						}

					};

					muMagThread.start();


				}
			}
		});

	}



	public static void main(String[] args) {
		new Main();
	}

	public void writeParameters(){

		//================================== writing the parameters and settings to the file:

		File f = new File(rootDirectory);

		if(!f.exists())
			f.mkdir();
	
		try{ PrintWriter pw = new  PrintWriter(new  FileWriter(parametersFile));

		pw.println("Parameters and settings:");
		pw.println();
		pw.println("A: "+model.A+" J/m");
		pw.println("Ms: "+model.Ms+" A/m");
		if(model.CubicAni)
			pw.println("Cubic Anisotropy: (x,y, z axes)");
		else if(model.uniAniAxis.el[0]==1)
			pw.println("Uniaxial Anisotropy: x axis)");
		else if(model.uniAniAxis.el[1]==1)
			pw.println("Uniaxial Anisotropy: y axis)");
		else if(model.uniAniAxis.el[2]==1)
			pw.println("Uniaxial Anisotropy: z axis)");

		pw.println("K1: "+model.K1+" J/m3");
		pw.println("K2: "+model.K2+" J/m3");

		pw.println("Lx: "+model.Lx+" nm");
		pw.println("Ly: "+model.Ly+" nm");
		pw.println("Lz: "+model.Lz+" nm");
		pw.println("Nx: "+model.Nx);
		pw.println("Ny "+model.Ny);
		pw.println("Nz: "+model.Nz);
		pw.println("dx: "+model.dxnm+" nm");
		pw.println("dy "+model.dynm+" nm");
		pw.println("dz: "+model.dznm+" nm");
		if(model.LLG){
			pw.println("Dynamic: LLG");
			pw.println("      gamma: "+model.gamma);
			pw.println("      alpha: "+model.alpha);
			pw.println("      dt: "+model.dt+" picosecons");
		}
		else{
			pw.println("Dynamic: Rotation");
			pw.println("      betta: "+model.betta);
		}
		if(model.PreviousInitial)
			pw.println("Previous result from file "+rootDirectory+prevmagFile+" was used as initial values.");
		else{
			pw.println("New distribution with the following parameters was used as initial values:");
			pw.println("      Mean of polar angle for m (Theta): "+thetaMean);
			pw.println("      Random range of polar angle: "+thetaRandomRange);
			pw.println("      Mean of azimuthal angle for m ( Phi): "+phiMean);
			pw.println("      Random range of azimuthal angle: "+phiRandomRange);
		}
		if(model.ModifEex) 
			{pw.println("Exchange Energt modified:");
		pw.println("Critical angle for exchange energy modification: "+thetaCx+" degree");}
		if(model.ModifEd)
			{pw.println("Exchange Energt modified:");
		pw.println("Critical angle for demagnetizing energy modification: "+thetaCd+" degree");}
		pw.println("Maximum Error: "+errorMax+" degree");


		pw.close();
		}

		catch (IOException e) { } 
		
		//================================== writing the parameters and settings to tGUI text area:
		out.println("============================================");
		out.println(" Run number "+runCounter);
		out.println("============================================");
		
		out.println(" Parameters and settings:");
		out.println();
		out.println(" Lx: "+model.Lx+" nm");
		out.println(" Ly: "+model.Ly+" nm");
		out.println(" Lz: "+model.Lz+" nm");
		out.println(" Nx: "+model.Nx);
		out.println(" Ny "+model.Ny);
		out.println(" Nz: "+model.Nz);
		out.println(" dx: "+model.dxnm+" nm");
		out.println(" dy "+model.dynm+" nm");
		out.println(" dz: "+model.dznm+" nm");
		out.println(" A: "+model.A+" J/m");
		out.println(" Ms: "+model.Ms+" A/m");
		if(model.CubicAni)
			out.println(" Cubic Anisotropy: (x,y, z axes)");
		else if(model.uniAniAxis.el[0]==1)
			out.println(" Uniaxial Anisotropy: x axis)");
		else if(model.uniAniAxis.el[1]==1)
			out.println(" Uniaxial Anisotropy: y axis)");
		else if(model.uniAniAxis.el[2]==1)
			out.println(" Uniaxial Anisotropy: z axis)");

		out.println(" K1: "+model.K1+" J/m3");
		out.println(" K2: "+model.K2+" J/m3");

		
		if(model.LLG){
			out.println(" Dynamic: LLG");
			out.println("      gamma: "+model.gamma);
			out.println("      alpha: "+model.alpha);
			out.println("      dt: "+model.dt+" picosecons");
		}
		else{
			out.println(" Dynamic: Rotation:");
			out.println("      betta: "+model.betta);
		}
		if(model.PreviousInitial)
			out.println(" Previous result from file "+rootDirectory+prevmagFile+" was used as initial values.");
		else{
			out.println(" New distribution with the following parameters was used as initial values:");
			out.println("      Mean of polar angle for m (Theta): "+thetaMean);
			out.println("      Random range of polar angle: "+thetaRandomRange);
			out.println("      Mean of azimuthal angle for m ( Phi): "+phiMean);
			out.println("      Random range of azimuthal angle: "+phiRandomRange);
		}
		if(model.ModifEex) 
			{out.println(" Exchange Energt modified:");
		out.println(" Critical angle for exchange energy modification: "+thetaCx+" degree");}
		
		if(model.ModifEd) 
			{out.println(" Exchange Energt modified:");
		out.println(" Critical angle for demagnetizing energy modification: "+thetaCd+" degree");
			}
		out.println(" Maximum Error: "+errorMax+" degree");

	}




}

