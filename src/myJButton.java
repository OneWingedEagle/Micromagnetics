import java.awt.Font;
import javax.swing.JButton;


public class myJButton extends JButton {
	public Font tfFont = new Font("Times New Roman", 1, 14);
	myJButton(){
	super();	
	
	setFont(tfFont);
	}
	public myJButton(String st, int alignment ){
		super();	
		setText(st);
		setHorizontalAlignment(alignment);
		setFont(tfFont);
		}
	public myJButton(String str){
		super();	
		setText(str);
		setFont(tfFont);
		}
}
