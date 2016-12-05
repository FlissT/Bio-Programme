import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

w = Gtk.Window()
data = Gtk.ListStore(int, str)
data.append([1, "test"])

combo = Gtk.ComboBox.new_with_model_and_entry(data)
renderer = Gtk.CellRendererText()
combo.pack_start(renderer, True)
combo.add_attribute(renderer, "text", 0)
combo.set_model(data)
combo.set_entry_text_column(1)
w.add(combo)

def on_file_selected(entry):
    file = open(self.fbox.get_filename())
    data.append([file])

w.connect("delete-event", Gtk.main_quit)


w.show_all()

Gtk.main()

