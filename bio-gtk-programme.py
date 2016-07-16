import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC


class LnFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file("seqlen-page-glade.glade")
        self.sl_box = self.builder.get_object("SeqLen-box")
        self.add(self.sl_box)
        self.show_all()
    
        self.entry = self.builder.get_object("SeqLen-entry")
        self.button = self.builder.get_object("SeqLen-button")
        self.label = self.builder.get_object("SeqLen-result")
        self.button.connect("clicked", self.clicked_callback)
        self.entry.connect("activate", self.enter_callback)
        

    def seq_len(self): #gives the length of the sequence
        self.sl_entry = self.builder.get_object("SeqLen-entry")
        entry_text = Seq(self.sl_entry.get_text(), IUPAC.unambiguous_dna)
        len_result = len(entry_text)
        #self.board_label.set_text("Result")
        self.label.set_text(str(len_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_len()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.seq_len()
        

class CmFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.builder = Gtk.Builder()
        self.builder.add_from_file("seqcomp-page-glade.glade")
        self.cm_box = self.builder.get_object("SeqComp-box")
        self.add(self.cm_box)
        self.show_all()

        self.entry1 = self.builder.get_object("SeqComp-entry1")
        self.entry2 = self.builder.get_object("SeqComp-entry2")
        self.button = self.builder.get_object("SeqComp-button")
        self.label = self.builder.get_object("SeqComp-result")
        self.button.connect("clicked", self.clicked_callback)
        
        self.entry1.connect("activate", self.enter_callback)
        self.entry2.connect("activate", self.enter_callback)


    def seq_comp(self): #compares two sequences
        self.cm_entry1 = self.builder.get_object("SeqComp-entry1")
        self.cm_entry2 = self.builder.get_object("SeqComp-entry2")
        entry_text1 = Seq(self.cm_entry1.get_text(), IUPAC.unambiguous_dna)
        entry_text2 = Seq(self.cm_entry2.get_text(), IUPAC.unambiguous_dna)
        comp_result = ((entry_text1) == (entry_text2))
        #self.board_label.set_text("Result")
        self.label.set_text(str(comp_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.seq_comp()
        

class GcFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.builder = Gtk.Builder()
        self.builder.add_from_file("gc-page-glade.glade")
        self.gc_box = self.builder.get_object("GC-box")
        self.add(self.gc_box)
        self.show_all()

        self.entry = self.builder.get_object("GC-entry")
        self.button = self.builder.get_object("GC-button")
        self.label = self.builder.get_object("GC-result")
        self.button.connect("clicked", self.clicked_callback)
        self.entry.connect("activate", self.enter_callback)
        

    def gc_content(self): #calculates gc % of sequence
        self.gc_entry = self.builder.get_object("GC-entry")
        entry_text = Seq(self.gc_entry.get_text(), IUPAC.unambiguous_dna)
        gc_result = GC(entry_text)
        #self.board_label.set_text("Result")
        self.label.set_text(str(gc_result))
        print(str(gc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.gc_content()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.gc_content()
        
class RcFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file("revcomp-page-glade.glade")
        self.rc_box = self.builder.get_object("RevComp-box")
        self.add(self.rc_box)
        self.show_all()

        self.entry = self.builder.get_object("RevComp-entry")
        self.button = self.builder.get_object("RevComp-button")
        self.label = self.builder.get_object("RevComp-result")
        self.button.connect("clicked", self.clicked_callback)
        self.entry.connect("activate", self.enter_callback)
       

    def rev_comp(self): #gives the reverse complement of DNA sequence
        self.rc_entry = self.builder.get_object("RevComp-entry")
        entry_text = Seq(self.rc_entry.get_text(), IUPAC.unambiguous_dna)
        rc_result = entry_text.reverse_complement()
        #self.board_label.set_text("Result")
        self.label.set_text(str(rc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.rev_comp()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.rev_comp()


class TrFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.builder = Gtk.Builder()
        self.builder.add_from_file("transl-page-glade.glade")
        self.tr_box = self.builder.get_object("Transl-box")
        self.add(self.tr_box)
        self.show_all()

        self.entry = self.builder.get_object("Transl-entry")
        self.button = self.builder.get_object("Transl-button")
        self.label = self.builder.get_object("Transl-result")
        self.button.connect("clicked", self.clicked_callback)
        self.entry.connect("activate", self.enter_callback)
        

    def translation(self): #translates sequence into protein
        self.tr_entry = self.builder.get_object("Transl-entry")
        entry_text = Seq(self.tr_entry.get_text(), IUPAC.unambiguous_dna)
        mrna = entry_text.transcribe()
        tr_result = mrna.translate()
        #self.board_label.set_text("Result")
        self.label.set_text(str(tr_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.translation()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.translation()


class EnFrame(Gtk.Bin):
    def __init__(self):
        Gtk.Bin.__init__(self)

        self.builder = Gtk.Builder()
        self.builder.add_from_file("entrez-page-glade.glade")
        self.en_box = self.builder.get_object("Entrez-box")
        self.add(self.en_box)
        self.show_all()

        self.entry = self.builder.get_object("Entrez-entry")
        self.button = self.builder.get_object("Entrez-button")
        self.label = self.builder.get_object("Entrez-result")
        self.button.connect("clicked", self.clicked_callback)
        self.fbox = self.builder.get_object("Entrez-file")
        #self.fbox.connect("file_set", self.on_file_selected)        
       
        self.entry.connect("activate", self.enter_callback)
        

    #def on_file_selected(self, entry):
        #print("You chose ", self.entry.get_filename())
        
    def entrez_db(self): #finds details from Entrez database
        self.tr_entry = self.builder.get_object("Entrez-entry")
        entry_text = Entrez.efetch(db = "nucleotide", id = [self.entry.get_text()], rettype = "fasta")
        en_result = SeqIO.read(entry_text, "fasta")
        #self.board_label.set_text("Result")
        self.label.set_text(str("Written to file 'bio-gtk-entrez.txt'"))
        SeqIO.write(en_result, "bio-gtk-entrez.txt", "fasta")
        print(en_result.id)
        print(en_result.seq)
        
    def clicked_callback(self, button): #runs function when button is clicked
        self.entrez_db()
        
    def enter_callback(self, widget): #runs function when enter pressed
        self.entrez_db()


class MyWindow(Gtk.Window):
    
    def __init__(self):
        Gtk.Window.__init__(self)
        self.set_size_request(300, 100)
        
        self.set_title("Bio Programme")
        self.msg_label = Gtk.Label()
        self.board_label = Gtk.Label("Bio")
        self.entry = Gtk.Entry()

        vbox = Gtk.VBox()
        stack = Gtk.Stack()
        stack.set_transition_type(Gtk.StackTransitionType.SLIDE_LEFT_RIGHT)
        stack.set_transition_duration(250)
        switcher = Gtk.StackSwitcher()
        switcher.set_stack(stack)
        vbox.add(switcher)
        vbox.add(stack)

        sl_box = LnFrame()
        cm_box = CmFrame()
        gc_box = GcFrame()
        rc_box = RcFrame()
        tr_box = TrFrame()
        en_box = EnFrame()

        stack.add_titled(sl_box, "sl", "Sequence Length")
        stack.add_titled(cm_box, "cm", "Sequence Comparison")
        stack.add_titled(gc_box, "gc", "Calculate GC")
        stack.add_titled(rc_box, "rc", "Reverse Complement")
        stack.add_titled(tr_box, "tr", "Translate Sequence")
        stack.add_titled(en_box, "en", "Entrez")
        self.add(vbox)
        
        self.show_all()
        self.connect("delete-event", self.on_quit)
        

    def on_quit(self, widget, event):
        Gtk.main_quit()

       
window = MyWindow()

Gtk.main()
