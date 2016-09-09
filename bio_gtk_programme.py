#CRASHES WHEN TRYING TO GET LARGE FILES FROM ENTREZ e.g. whole chromosomes, fine with single genes

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from Bio import Entrez, SeqIO, Medline
Entrez.email = "A.N.Other@example.com" 
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC


class LnFrame(Gtk.Bin):
    def __init__(self, open_sequences, seq_liststore):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        self.seq_liststore = seq_liststore
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file("seqlen-page-glade.glade")
        self.sl_box = self.builder.get_object("SeqLen-box")
        self.add(self.sl_box)
    
        self.cbox = self.builder.get_object("SeqLen-cbox")
        self.button = self.builder.get_object("SeqLen-button")
        self.label = self.builder.get_object("SeqLen-result")
        self.button.connect("clicked", self.clicked_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(seq_liststore)
        
    def seq_len(self): #gives the length of the sequence
        iterator = self.cbox.get_active_iter()
        seq_id = self.seq_liststore.get_value(iterator, 0)
        seq = (self.open_sequences[seq_id])
        len_result = len(seq)
        self.label.set_text(str(len_result))
        
    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_len()
        

class CmFrame(Gtk.Bin):
    def __init__(self, open_sequences, seq_liststore):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        self.seq_liststore = seq_liststore

        self.builder = Gtk.Builder()
        self.builder.add_from_file("seqcomp-page-glade.glade")
        self.cm_box = self.builder.get_object("SeqComp-box")
        self.add(self.cm_box)

        self.cbox1 = self.builder.get_object("SeqComp-cbox1")
        self.cbox2 = self.builder.get_object("SeqComp-cbox2")
        self.button = self.builder.get_object("SeqComp-button")
        self.label = self.builder.get_object("SeqComp-result")
        self.button.connect("clicked", self.clicked_callback)
        
        renderer = Gtk.CellRendererText()
        self.cbox1.pack_start(renderer, True)
        self.cbox1.add_attribute(renderer, "text", 0)
        self.cbox1.set_model(seq_liststore)
        self.cbox2.pack_start(renderer, True)
        self.cbox2.add_attribute(renderer, "text", 0)
        self.cbox2.set_model(seq_liststore)
        
    def seq_comp(self): #compares two sequences
        iterator1 = self.cbox1.get_active_iter()
        seq_id1 = self.seq_liststore.get_value(iterator1, 0)
        seq1 = self.open_sequences[seq_id1]
        iterator2 = self.cbox2.get_active_iter()
        seq_id2 = self.seq_liststore.get_value(iterator2, 0)
        seq2 = self.open_sequences[seq_id2]
        comp_result = ((seq1) == (seq2))
        self.label.set_text(str(comp_result))


    def clicked_callback(self, button): #runs function when button is clicked
        self.seq_comp()
        

class GcFrame(Gtk.Bin):
    def __init__(self, open_sequences, seq_liststore):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        self.seq_liststore = seq_liststore

        self.builder = Gtk.Builder()
        self.builder.add_from_file("gc-page-glade.glade")
        self.gc_box = self.builder.get_object("GC-box")
        self.add(self.gc_box)

        self.cbox = self.builder.get_object("GC-cbox")
        self.button = self.builder.get_object("GC-button")
        self.label = self.builder.get_object("GC-result")
        self.button.connect("clicked", self.clicked_callback)
        
        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(seq_liststore)
        
    def gc_content(self): #calculates gc % of sequence
        iterator = self.cbox.get_active_iter()
        seq_id = self.seq_liststore.get_value(iterator, 0)
        seq = self.open_sequences[seq_id]
        gc_result = round(GC(seq), 3)
        self.label.set_text(str(gc_result))

    def clicked_callback(self, button): #runs function when button is clicked
        self.gc_content()
        

class RcFrame(Gtk.Bin):
    def __init__(self, open_sequences, seq_liststore):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        self.seq_liststore = seq_liststore
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file("revcomp-page-glade.glade")
        self.rc_box = self.builder.get_object("RevComp-box")
        self.add(self.rc_box)

        self.cbox = self.builder.get_object("RevComp-cbox")
        self.button = self.builder.get_object("RevComp-button")
        self.label = self.builder.get_object("RevComp-result")
        self.button.connect("clicked", self.clicked_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(seq_liststore)

    
    def rev_comp(self): #gives the reverse complement of DNA sequence
        iterator = self.cbox.get_active_iter()
        seq_id = self.seq_liststore.get_value(iterator, 0)
        seq = self.open_sequences[seq_id]
        rc_result = seq.reverse_complement()                
        
        tv = Gtk.TextView()
        tv.get_buffer().set_text(str(rc_result))
        tv.set_editable(False)          
        sw = Gtk.ScrolledWindow()
        sw.set_size_request(300,200)
        sw.add(tv)
        w = Gtk.Window()                                                                                                                                                                                                                                                
        w.add(sw)
        w.show_all()

    def clicked_callback(self, button): #runs function when button is clicked
        self.rev_comp()


class TrFrame(Gtk.Bin):
    def __init__(self, open_sequences, seq_liststore):
        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        self.seq_liststore = seq_liststore

        self.builder = Gtk.Builder()
        self.builder.add_from_file("transl-page-glade.glade")
        self.tr_box = self.builder.get_object("Transl-box")
        self.add(self.tr_box)

        self.cbox = self.builder.get_object("Transl-cbox")
        self.button1 = self.builder.get_object("Transl-button")
        self.button2 = self.builder.get_object("Codon-button")
        self.label = self.builder.get_object("Transl-result")
        self.button1.connect("clicked", self.clicked_callback1)
        self.button2.connect("clicked", self.clicked_callback2)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(seq_liststore)
        
    def translation(self): #translates sequence into protein
        iterator = self.cbox.get_active_iter()
        seq_id = self.seq_liststore.get_value(iterator, 0)
        seq = self.open_sequences[seq_id]
        mrna = seq.transcribe()
        tr_result = mrna.translate()

        if len(mrna) % 3 != 0:
            self.label.set_text("Warning, incomplete codon")

        tv = Gtk.TextView()
        tv.get_buffer().set_text(str(tr_result))
        tv.set_editable(False)          
        sw = Gtk.ScrolledWindow()
        sw.set_size_request(200,100)
        sw.add(tv)
        w = Gtk.Window()
        w.add(sw)
        w.show_all()

    def clicked_callback1(self, button1): #runs function when button is clicked
        self.translation()

    def codon(self):
        from Bio.Data import CodonTable
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.label.set_text(str(standard_table))

    def clicked_callback2(self, button2):
        self.codon()
        

class OsFrame(Gtk.Bin): #opens sequences for later use
    def __init__(self, open_sequences, seq_liststore):

        Gtk.Bin.__init__(self)

        self.open_sequences = open_sequences
        self.seq_liststore = seq_liststore

        self.builder = Gtk.Builder()
        self.builder.add_from_file("openseq-page.glade")
        self.os_box = self.builder.get_object("Open-box")
        self.add(self.os_box)

        self.cbox = self.builder.get_object("Open-cbox")
        self.fbox = self.builder.get_object("Open-file")
        self.fbox.connect("file_set", self.on_file_selected)

        self.entry = self.builder.get_object("Entrez-entry")
        self.button = self.builder.get_object("Entrez-button")
        self.label_e = self.builder.get_object("Entrez-label")
        self.label = self.builder.get_object("Entrez-result")
        self.label1 = self.builder.get_object("Entrez-result1")
        self.button.connect("clicked", self.clicked_callback)

        renderer = Gtk.CellRendererText()
        self.cbox.pack_start(renderer, True)
        self.cbox.add_attribute(renderer, "text", 0)
        self.cbox.set_model(seq_liststore)
        self.cbox.set_active(0)
                
    def on_file_selected(self, entry): #opens a file and adds it to list store
        file = open(self.fbox.get_filename())
        for seq_record in SeqIO.parse(file, "fasta"):
            self.open_sequences[seq_record.id] = seq_record.seq # Add the id and sequence to the sequences dict
            self.seq_liststore.append([seq_record.id]) # Add the id to the liststore so we can look it up later


    def entrez_db(self): #finds details from Entrez database
        entry_text = Entrez.efetch(db = "nucleotide", id = [self.entry.get_text()], rettype = "fasta")
        en_result = SeqIO.read(entry_text, "fasta")
        self.label_e.set_text(str("Written to file (Entrez ID)"))
        SeqIO.write(en_result, en_result.id, "fasta")
        self.label.set_text(str("Open file"))
        self.label1.set_text(str(en_result.description))


    def clicked_callback(self, button): #runs function when button is clicked
        self.entrez_db()

class PbFrame(Gtk.Bin): #opens sequences for later use
    def __init__(self):

        Gtk.Bin.__init__(self)

        self.builder = Gtk.Builder()
        self.builder.add_from_file("pubmed-page-glade.glade")
        self.pb_box = self.builder.get_object("Pbmd-box")
        self.add(self.pb_box)

        self.entry = self.builder.get_object("Pbmd-entry")
        self.button = self.builder.get_object("Pbmd-button")
        self.label1 = self.builder.get_object("Pbmd-label")
        self.label = self.builder.get_object("Pbmd-result")
        self.button.connect("clicked", self.clicked_callback)
                
    def pbmd_search(self): #searches pubmed database
        handle = Entrez.egquery(term=self.entry.get_text())
        record = Entrez.read(handle)
        for row in record["eGQueryResult"]:
            if row["DbName"]=="pubmed":
                print(row["Count"])

        handle = Entrez.esearch(db="pubmed", term=self.entry.get_text(), retmax=1000)
        record = Entrez.read(handle)
        idlist = record["IdList"]
        
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        records = list(records)
        
        for record in records: #want to get this bit in a Gtk.TextView
            print("title:", record.get("TI", "?"))
            print("authors:", record.get("AU", "?"))
            print("source:", record.get("SO", "?"))
            print("")
            

    def clicked_callback(self, button): #runs function when button is clicked
        self.pbmd_search()
        
        
class MyWindow(Gtk.Window):
    
    def __init__(self):
        Gtk.Window.__init__(self)
        self.set_size_request(300, 300)
        
        vbox = Gtk.VBox()       
        nbook = Gtk.Notebook()
        nbook.set_tab_pos(Gtk.PositionType.LEFT)
        nbook.set_scrollable(True)
        vbox.add(nbook)

        self.open_sequences = {}
        self.seq_liststore = Gtk.ListStore(str)

        self.set_title("Bio Programme")
                          
        os_box = OsFrame(self.open_sequences, self.seq_liststore)
        pb_box = PbFrame()
        sl_box = LnFrame(self.open_sequences, self.seq_liststore)
        cm_box = CmFrame(self.open_sequences, self.seq_liststore)
        gc_box = GcFrame(self.open_sequences, self.seq_liststore)
        rc_box = RcFrame(self.open_sequences, self.seq_liststore)
        tr_box = TrFrame(self.open_sequences, self.seq_liststore)
        
        nbook.append_page(os_box)
        nbook.append_page(pb_box)
        nbook.append_page(sl_box)
        nbook.append_page(cm_box)
        nbook.append_page(gc_box)
        nbook.append_page(rc_box)
        nbook.append_page(tr_box)
        
        nbook.set_tab_label_text(os_box, "Get Sequences")
        nbook.set_tab_label_text(pb_box, "Search Pubmed")
        nbook.set_tab_label_text(sl_box, "Sequence Length")
        nbook.set_tab_label_text(cm_box, "Sequence Comparison")
        nbook.set_tab_label_text(gc_box, "GC Content")
        nbook.set_tab_label_text(rc_box, "Reverse Complement")
        nbook.set_tab_label_text(tr_box, "Translate Sequence")
        self.add(vbox)
        
        self.show_all()
        self.connect("delete-event", self.on_quit)
        

    def on_quit(self, widget, event):
        Gtk.main_quit()

       
window = MyWindow()

Gtk.main()
