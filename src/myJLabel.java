import java.awt.Font;

import javax.swing.JLabel;
public class myJLabel extends JLabel {
	Font tfFont = new Font("Arial", 0, 13);
	myJLabel(){
	super();	
	
	setFont(tfFont);
	}
	public myJLabel(String st, int alignment ){
		super();	
		setText(st);
		setHorizontalAlignment(alignment);
		setFont(tfFont);
		}
	
	public myJLabel(String str){
		super();	
		setText(str);
		setFont(tfFont);
		}
}
