import java.awt.Dimension;
import java.awt.Font;
import javax.swing.JTextField;

public class myJTextField extends JTextField {
	Font tfFont = new Font("Arial", 0, 14);
	myJTextField(){
	super();	
	
	setFont(tfFont);
	}
	public myJTextField(int size ){
		super();	
		setSize(new Dimension(10,size));
		setFont(tfFont);
		}
	public myJTextField(String str ){
		super();	
		setText(str);
		setFont(tfFont);
		}
	public myJTextField(int size, String str ){
		super();	
		setText(str);
		setSize(new Dimension(10,size));
		setFont(tfFont);
		}
}
