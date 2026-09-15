#include "confedit.hpp"

// The keywords, types, documentation and snippets used to be hard-coded here,
// in four init_*() functions of some 680 lines. They now live in rockable.lang,
// which the terminal editor of prepro/rockedit reads too: documenting a new
// keyword is a matter of editing a text file, and the two editors cannot drift
// apart.
rockable_lang::Language language;

void init_language(const char* argv0) {
  const std::string file = rockable_lang::findFile(argv0);
  if (file.empty() || !language.load(file)) {
    fl_alert(
        "rockable.lang was not found.\n\n"
        "Syntax highlighting, documentation and snippets are therefore disabled.\n\n"
        "Point the ROCKABLE_LANG environment variable at it, or copy it next to\n"
        "the confedit executable.");
  }
}

// Returns the stream's current offset, clamped to [0, length]. std::istream
// returns -1 from tellg() once EOF has been reached (e.g. after extracting the
// last token of a buffer without a trailing separator); without this clamp the
// computed offsets go negative and we write style[] out of bounds.
static int stream_pos(std::istream& is, int length) {
  std::streampos p = is.tellg();
  if (p == std::streampos(-1)) return length;
  int pos = static_cast<int>(p);
  if (pos < 0) return 0;
  if (pos > length) return length;
  return pos;
}

void style_parse(const char* text, char* style, int length) {
  for (int i = 0; i < length; i++) {
    if (text[i] != '\n') {
      style[i] = 'A';
    }
  }
  std::istringstream is(std::string(text, length));

  std::string token;
  while (is >> token) {
    int to = stream_pos(is, length);
    int from = to - static_cast<int>(token.size());
    if (from < 0) from = 0;

    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      getline(is, token);  // ignore the rest of the current line
      int lineEnd = stream_pos(is, length);
      for (int pos = from; pos < lineEnd; pos++) {
        style[pos] = 'B';
      }
      continue;
    }

    const rockable_lang::Kind kind = language.kindOf(token);
    if (kind == rockable_lang::Kind::Keyword) {
      for (int pos = from; pos < to; pos++) {
        style[pos] = 'D';
      }
    } else if (kind == rockable_lang::Kind::Type) {
      for (int pos = from; pos < to; pos++) {
        style[pos] = 'C';
      }
    }
  }  // while
}

void style_init(void) {
  char* style = new char[textbuf->length() + 1];
  char* text = textbuf->text();

  memset(style, 'A', textbuf->length());
  style[textbuf->length()] = '\0';

  if (!stylebuf) {
    stylebuf = new Fl_Text_Buffer(textbuf->length());
  }

  style_parse(text, style, textbuf->length());

  stylebuf->text(style);
  delete[] style;
  free(text);
}

void style_unfinished_cb(int, void*) {}

void style_update(int pos,                      // I - Position of update
                  int nInserted,                // I - Number of inserted chars
                  int nDeleted,                 // I - Number of deleted chars
                  int /*nRestyled*/,            // I - Number of restyled chars
                  const char* /*deletedText*/,  // I - Text that was deleted
                  void* cbArg) {                // I - Callback data
  int start,                                    // Start of text
      end;                                      // End of text
  char *style,                                  // Style data
      *text;                                    // Text data

  // If this is just a selection change, just unselect the style buffer...
  if (nInserted == 0 && nDeleted == 0) {
    stylebuf->unselect();
    return;
  }

  // Track changes in the text buffer...
  if (nInserted > 0) {
    // Insert characters into the style buffer...
    style = new char[nInserted + 1];
    memset(style, 'A', nInserted);
    style[nInserted] = '\0';

    stylebuf->replace(pos, pos + nDeleted, style);
    delete[] style;
  } else {
    // Just delete characters in the style buffer...
    stylebuf->remove(pos, pos + nDeleted);
  }

  // Select the area that was just updated to avoid unnecessary
  // callbacks...
  stylebuf->select(pos, pos + nInserted - nDeleted);

  // Re-parse only the line(s) touched by the change. The Rockable conf format
  // has no multi-line constructs (comments end at the line), so a line's style
  // depends solely on its own content; there is no need to reparse the rest of
  // the buffer, which would be very slow on large files.
  start = textbuf->line_start(pos);
  end = textbuf->line_end(pos + nInserted);
  text = textbuf->text_range(start, end);
  style = stylebuf->text_range(start, end);

  style_parse(text, style, end - start);

  stylebuf->replace(start, end, style);
  ((Fl_Text_Editor*)cbArg)->redisplay_range(start, end);

  free(text);
  free(style);
}

EditorWindow::EditorWindow(int w, int h, const char* t) : Fl_Double_Window(w, h, t) {

  // replace dialog
  replace_dlg = new Fl_Window(300, 105, "Replace");
  replace_find = new Fl_Input(80, 10, 210, 25, "Find:");
  replace_find->align(FL_ALIGN_LEFT);

  replace_with = new Fl_Input(80, 40, 210, 25, "Replace:");
  replace_with->align(FL_ALIGN_LEFT);

  replace_all = new Fl_Button(10, 70, 90, 25, "Replace All");
  replace_all->callback((Fl_Callback*)replall_cb, this);

  replace_next = new Fl_Return_Button(105, 70, 120, 25, "Replace Next");
  replace_next->callback((Fl_Callback*)replace2_cb, this);

  replace_cancel = new Fl_Button(230, 70, 60, 25, "Cancel");
  replace_cancel->callback((Fl_Callback*)replcan_cb, this);
  replace_dlg->end();
  replace_dlg->set_non_modal();

  // keyword dialog
  keyword_dlg = new Fl_Window(300, 400, "Keywords");

  keyword_browser = new Fl_Browser(10, 10, 280, 380 - 25 - 5);
  keyword_browser->type(FL_HOLD_BROWSER);

  for (const rockable_lang::Snippet& snippet : language.snippets()) {
    keyword_browser->add(snippet.title.c_str());
  }

  keyword_insert = new Fl_Button(230 - 60 - 10, 365, 60, 25, "Insert");
  keyword_insert->callback((Fl_Callback*)keyword_insert_cb, this);

  keyword_close = new Fl_Button(230, 365, 60, 25, "Close");
  keyword_close->callback((Fl_Callback*)keyword_close_cb, this);

  keyword_dlg->end();
  keyword_dlg->set_non_modal();

  // add particle dialog
  add_particle_dlg = new addParticleDialog(400, 150, "Add a particle");

  //  ========
  editor = 0;
  *search = (char)0;
  wrap_mode = 0;
  line_numbers = 0;
}

EditorWindow::~EditorWindow() {
  delete replace_dlg;
  delete keyword_dlg;
  delete add_particle_dlg;
}

void doc_selection_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Buffer* buf = e->editor->buffer();
  if (buf->selected()) {
    const std::string* doc = language.doc(buf->selection_text());
    if (doc != nullptr) {
      fl_message("%s", doc->c_str());
    }
  }
}

void deg_to_rad_cb(Fl_Widget*, void* v) {
  const char* deg = fl_input("Value in degree: ", "0");
  char radtxt[100];
  snprintf(radtxt, 100, "%f", atof(deg) * M_PI / 180.0);

  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Buffer* buf = e->editor->buffer();
  buf->insert(e->editor->insert_position(), radtxt);
}

void quat_axis_angle_cb(Fl_Widget*, void* v) {
  const char* input = fl_input("Enter axis vector followed by an angle (degree): ", "1 0 0 90");
  std::stringstream ss(input);
  double axis_x, axis_y, axis_z, angleDeg;
  ss >> axis_x >> axis_y >> axis_z >> angleDeg;
  double s, vx, vy, vz;
  double half_angle = 0.5 * angleDeg * M_PI / 180.;
  double len = sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
  s = cos(half_angle);
  vx = sin(half_angle) * axis_x / len;
  vy = sin(half_angle) * axis_y / len;
  vz = sin(half_angle) * axis_z / len;
  char txt[100];
  snprintf(txt, 100, "%f  %f %f %f", s, vx, vy, vz);

  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Buffer* buf = e->editor->buffer();
  buf->insert(e->editor->insert_position(), txt);
}

void add_particle_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->add_particle_dlg->show();
}

int check_save(void) {
  if (!changed) return 1;

  int r = fl_choice(
      "The current file has not been saved.\n"
      "Would you like to save it now?",
      "Cancel", "Save", "Don't Save");

  if (r == 1) {
    save_cb();  // Save the file...
    return !changed;
  }

  return (r == 2) ? 1 : 0;
}

void load_file(const char* newfile, int ipos) {
  loading = 1;
  int insert = (ipos != -1);
  changed = insert;
  if (!insert) {
    filename[0] = '\0';
  }
  int r;
  if (!insert) {
    r = textbuf->loadfile(newfile);
  } else {
    r = textbuf->insertfile(newfile, ipos);
  }
  changed = changed || textbuf->input_file_was_transcoded;
  if (r) {
    fl_alert("Error reading from file \'%s\':\n%s.", newfile, strerror(errno));
  } else if (!insert) {
    snprintf(filename, sizeof(filename), "%s", newfile);
  }
  loading = 0;
  textbuf->call_modify_callbacks();
}

void save_file(const char* newfile) {
  if (textbuf->savefile(newfile)) {
    fl_alert("Error writing to file \'%s\':\n%s.", newfile, strerror(errno));
  } else {
    snprintf(filename, sizeof(filename), "%s", newfile);
  }
  changed = 0;
  textbuf->call_modify_callbacks();
}

void copy_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Editor::kf_copy(0, e->editor);
}

void cut_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Editor::kf_cut(0, e->editor);
}

void delete_cb(Fl_Widget*, void*) { textbuf->remove_selection(); }

void linenumbers_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Menu_Bar* m = (Fl_Menu_Bar*)w;
  const Fl_Menu_Item* i = m->mvalue();
  if (i->value()) {
    e->editor->linenumber_width(line_num_width);  // enable
    e->editor->linenumber_size(e->editor->textsize());
  } else {
    e->editor->linenumber_width(0);  // disable
  }
  e->line_numbers = (i->value() ? 1 : 0);
  e->redraw();
}

void wordwrap_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Menu_Bar* m = (Fl_Menu_Bar*)w;
  const Fl_Menu_Item* i = m->mvalue();
  if (i->value()) {
    e->editor->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS, 0);
  } else {
    e->editor->wrap_mode(Fl_Text_Display::WRAP_NONE, 0);
  }
  e->wrap_mode = (i->value() ? 1 : 0);
  e->redraw();
}

void find_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* val;

  val = fl_input("Search String:", e->search);
  if (val != NULL) {
    // User entered a string - go find it!
    snprintf(e->search, sizeof(e->search), "%s", val);
    find2_cb(w, v);
  }
}

void find2_cb(Fl_Widget* w, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  if (e->search[0] == '\0') {
    // Search string is blank; get a new one...
    find_cb(w, v);
    return;
  }

  int pos = e->editor->insert_position();
  int found = textbuf->search_forward(pos, e->search, &pos);
  if (found) {
    // Found a match; select and update the position...
    textbuf->select(pos, pos + strlen(e->search));
    e->editor->insert_position(pos + strlen(e->search));
    e->editor->show_insert_position();
  } else {
    fl_alert("No occurrences of \'%s\' found!", e->search);
  }
}

void set_title(Fl_Window* w) {
  if (filename[0] == '\0') {
    snprintf(title, sizeof(title), "Untitled");
  } else {
    char* slash;
    slash = strrchr(filename, '/');
#ifdef WIN32
    if (slash == NULL) {
      slash = strrchr(filename, '\\');
    }
#endif
    if (slash != NULL) {
      snprintf(title, sizeof(title), "%s", slash + 1);
    } else {
      snprintf(title, sizeof(title), "%s", filename);
    }
  }

  if (changed) {
    size_t len = strlen(title);
    snprintf(title + len, sizeof(title) - len, " (modified)");
  }

  w->label(title);
}

void changed_cb(int, int nInserted, int nDeleted, int, const char*, void* v) {
  if ((nInserted || nDeleted) && !loading) {
    changed = 1;
  }
  EditorWindow* w = (EditorWindow*)v;
  set_title(w);
  update_status(w);
  if (loading) {
    w->editor->show_insert_position();
  }
}

void new_cb(Fl_Widget*, void*) {
  if (!check_save()) {
    return;
  }

  filename[0] = '\0';
  textbuf->select(0, textbuf->length());
  textbuf->remove_selection();
  changed = 0;
  textbuf->call_modify_callbacks();
}

void open_cb(Fl_Widget*, void*) {
  if (!check_save()) {
    return;
  }
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Open file");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  if (fnfc.show()) {
    return;
  }
  load_file(fnfc.filename(), -1);
}

void insert_cb(Fl_Widget*, void* v) {
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Insert file");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_FILE);
  if (fnfc.show()) {
    return;
  }
  EditorWindow* w = (EditorWindow*)v;
  load_file(fnfc.filename(), w->editor->insert_position());
}

void paste_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  Fl_Text_Editor::kf_paste(0, e->editor);
}

void close_cb(Fl_Widget*, void* v) {
  EditorWindow* w = (EditorWindow*)v;

  if (num_windows == 1) {
    if (!check_save()) {
      return;
    }
  }

  w->hide();
  w->editor->buffer(0);
  textbuf->remove_modify_callback(style_update, w->editor);
  textbuf->remove_modify_callback(changed_cb, w);
  Fl::delete_widget(w);

  num_windows--;
  if (!num_windows) {
    exit(0);
  }
}

void quit_cb(Fl_Widget*, void*) {
  if (changed && !check_save()) {
    return;
  }

  exit(0);
}

void replace_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->replace_dlg->show();
}

void keyword_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->keyword_dlg->show();
}

void replace2_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* find = e->replace_find->value();
  const char* replace = e->replace_with->value();

  if (find[0] == '\0') {
    // Search string is blank; get a new one...
    e->replace_dlg->show();
    return;
  }

  e->replace_dlg->hide();

  int pos = e->editor->insert_position();
  int found = textbuf->search_forward(pos, find, &pos);

  if (found) {
    // Found a match; update the position and replace text...
    textbuf->select(pos, pos + strlen(find));
    textbuf->remove_selection();
    textbuf->insert(pos, replace);
    textbuf->select(pos, pos + strlen(replace));
    e->editor->insert_position(pos + strlen(replace));
    e->editor->show_insert_position();
  } else {
    fl_alert("No occurrences of \'%s\' found!", find);
  }
}

void replall_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* find = e->replace_find->value();
  const char* replace = e->replace_with->value();

  find = e->replace_find->value();
  if (find[0] == '\0') {
    // Search string is blank; get a new one...
    e->replace_dlg->show();
    return;
  }

  e->replace_dlg->hide();

  e->editor->insert_position(0);
  int times = 0;

  // Loop through the whole string
  for (int found = 1; found;) {
    int pos = e->editor->insert_position();
    found = textbuf->search_forward(pos, find, &pos);

    if (found) {
      // Found a match; update the position and replace text...
      textbuf->select(pos, pos + strlen(find));
      textbuf->remove_selection();
      textbuf->insert(pos, replace);
      e->editor->insert_position(pos + strlen(replace));
      e->editor->show_insert_position();
      times++;
    }
  }

  if (times) {
    fl_message("Replaced %d occurrences.", times);
  } else {
    fl_alert("No occurrences of \'%s\' found!", find);
  }
}

void replcan_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->replace_dlg->hide();
}

void keyword_close_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->keyword_dlg->hide();
}

void keyword_insert_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  int line = e->keyword_browser->value();
  if (line > 0) {
    // The browser lists the snippets in the order of rockable.lang, so the
    // selected row is directly the snippet index.
    const size_t index = static_cast<size_t>(line - 1);
    if (index < language.snippets().size()) {
      const std::string toInsert = language.snippets()[index].body + "\n";
      textbuf->insert(e->editor->insert_position(), toInsert.c_str());
    }
  }
}

void save_cb() {
  if (filename[0] == '\0') {
    // No filename - get one!
    saveas_cb();
    return;
  } else {
    save_file(filename);
  }
}

void saveas_cb() {
  Fl_Native_File_Chooser fnfc;
  fnfc.title("Save File As?");
  fnfc.type(Fl_Native_File_Chooser::BROWSE_SAVE_FILE);
  if (fnfc.show()) {
    return;
  }
  save_file(fnfc.filename());
}

void view_cb(Fl_Widget*, void*) {
  Fl_Window* w = new_view();
  w->show();
}

void goto_start_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->editor->insert_position(0);
  e->editor->show_insert_position();
}

void goto_end_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  e->editor->insert_position(textbuf->length());
  e->editor->show_insert_position();
}

void goto_line_cb(Fl_Widget*, void* v) {
  EditorWindow* e = (EditorWindow*)v;
  const char* val = fl_input("Go to line:", "1");
  if (val == NULL) {
    return;
  }
  int line = atoi(val);
  if (line < 1) {
    line = 1;
  }
  // skip_lines(0, line - 1) returns the buffer position at the start of 'line'
  int pos = textbuf->skip_lines(0, line - 1);
  e->editor->insert_position(pos);
  e->editor->show_insert_position();
}

// ----- Status bar -----------------------------------------------------------

void update_status(EditorWindow* w) {
  if (!w || !w->status || !w->editor) {
    return;
  }
  int pos = w->editor->insert_position();
  int line = textbuf->count_lines(0, pos) + 1;
  int col = pos - textbuf->line_start(pos) + 1;
  char buf[128];
  snprintf(buf, sizeof(buf), "  Ln %d, Col %d  |  %d chars", line, col, textbuf->length());
  w->status->copy_label(buf);
}

int ConfTextEditor::handle(int e) {
  int r = Fl_Text_Editor::handle(e);
  update_status(win);
  return r;
}

// ----- Theme ----------------------------------------------------------------

int dark_mode = 1;  // dark is the default theme

void apply_theme_to_editor(Fl_Text_Editor* ed) {
  if (dark_mode) {
    ed->color(fl_rgb_color(30, 30, 30));
    ed->textcolor(fl_rgb_color(220, 220, 220));
    ed->selection_color(fl_rgb_color(60, 90, 140));
    ed->cursor_color(fl_rgb_color(80, 170, 255));
    ed->linenumber_bgcolor(fl_rgb_color(58, 58, 58));
    ed->linenumber_fgcolor(fl_rgb_color(140, 140, 140));
  } else {
    ed->color(FL_WHITE);
    ed->textcolor(fl_rgb_color(30, 30, 30));
    ed->selection_color(fl_rgb_color(180, 210, 255));
    ed->cursor_color(fl_rgb_color(0, 120, 215));
    ed->linenumber_bgcolor(fl_rgb_color(240, 240, 240));
    ed->linenumber_fgcolor(fl_rgb_color(150, 150, 150));
  }
}

// Theme the window "chrome" (menu bar + status bar) plus the editor itself.
void apply_theme_to_window(EditorWindow* e) {
  if (!e) {
    return;
  }
  Fl_Color bg, fg;
  if (dark_mode) {
    bg = fl_rgb_color(45, 45, 45);
    fg = fl_rgb_color(220, 220, 220);
  } else {
    bg = fl_rgb_color(238, 238, 238);
    fg = fl_rgb_color(40, 40, 40);
  }
  if (e->editor) {
    apply_theme_to_editor(e->editor);
  }
  if (e->menubar) {
    e->menubar->color(bg);
    e->menubar->textcolor(fg);
  }
  if (e->status) {
    e->status->color(bg);
    e->status->labelcolor(fg);
  }
  e->redraw();
}

void apply_theme() {
  if (dark_mode) {
    // Syntax colors tuned for a dark background
    styletable[0].color = fl_rgb_color(220, 220, 220);  // A - Plain
    styletable[1].color = fl_rgb_color(120, 190, 120);  // B - Comments
    styletable[2].color = fl_rgb_color(90, 200, 200);   // C - Types
    styletable[3].color = fl_rgb_color(110, 160, 255);  // D - Keywords
    Fl::background(45, 45, 45);
    Fl::background2(30, 30, 30);
    Fl::foreground(220, 220, 220);
  } else {
    styletable[0].color = FL_BLACK;       // A - Plain
    styletable[1].color = FL_DARK_GREEN;  // B - Comments
    styletable[2].color = FL_DARK_CYAN;   // C - Types
    styletable[3].color = FL_BLUE;        // D - Keywords
    Fl::background(238, 238, 238);
    Fl::background2(255, 255, 255);
    Fl::foreground(30, 30, 30);
  }
  Fl::reload_scheme();

  // Push the new look to every open editor view
  for (Fl_Window* win = Fl::first_window(); win; win = Fl::next_window(win)) {
    apply_theme_to_window(dynamic_cast<EditorWindow*>(win));
  }
}

void theme_light_cb(Fl_Widget*, void*) {
  dark_mode = 0;
  apply_theme();
}

void theme_dark_cb(Fl_Widget*, void*) {
  dark_mode = 1;
  apply_theme();
}

void ConfTextEditor::draw() {
  Fl_Text_Editor::draw();
  int lnw = linenumber_width();
  if (lnw <= 0) {
    return;  // gutter hidden -> no separator
  }
  // Vertical separator at the right edge of the line-number gutter
  fl_color(dark_mode ? fl_rgb_color(90, 90, 90) : fl_rgb_color(190, 190, 190));
  int sepx = x() + lnw - 1;
  fl_push_clip(x(), y(), w(), h());
  fl_line(sepx, y() + 1, sepx, y() + h() - 2);
  fl_pop_clip();
}

Fl_Window* new_view() {
  EditorWindow* w = new EditorWindow(with0, height0, title);

  const int menuH = 30;
  const int statusH = 24;

  w->begin();
  Fl_Menu_Bar* m = new Fl_Menu_Bar(0, 0, with0, menuH);
  m->copy(menuitems, w);
  m->box(FL_FLAT_BOX);  // honor color() instead of the scheme's gradient
  w->menubar = m;

  ConfTextEditor* ed = new ConfTextEditor(0, menuH, with0, height0 - menuH - statusH);
  ed->win = w;
  w->editor = ed;
  w->editor->textfont(TF);
  w->editor->textsize(TS);
  w->editor->cursor_style(Fl_Text_Display::SIMPLE_CURSOR);
  // w->editor->wrap_mode(Fl_Text_Editor::WRAP_AT_BOUNDS, 250);
  w->editor->buffer(textbuf);

  // Line numbers shown by default (matches the FL_MENU_VALUE flag on the menu toggle)
  w->editor->linenumber_width(line_num_width);
  w->editor->linenumber_size(TS);
  w->editor->linenumber_font(TF);
  w->editor->linenumber_align(FL_ALIGN_RIGHT);
  w->line_numbers = 1;

  w->editor->highlight_data(stylebuf, styletable, sizeof(styletable) / sizeof(styletable[0]), 'A', style_unfinished_cb,
                            0);
  w->plugDialogsWithEditor();

  // Status bar pinned to the bottom (editor is the resizable widget above it)
  w->status = new Fl_Box(0, height0 - statusH, with0, statusH, "Ln 1, Col 1");
  w->status->box(FL_FLAT_BOX);
  w->status->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
  w->status->labelsize(12);
  w->status->labelfont(FL_HELVETICA);

  apply_theme_to_window(w);  // colors for editor + menu bar + status (current theme)

  w->end();
  w->resizable(w->editor);
  w->size_range(300, 200);
  w->callback((Fl_Callback*)close_cb, w);

  textbuf->add_modify_callback(style_update, w->editor);
  textbuf->add_modify_callback(changed_cb, w);
  textbuf->call_modify_callbacks();
  num_windows++;
  return w;
}

void cb(const char* fname) { load_file(fname, -1); }

int main(int argc, char** argv) {
  Fl::scheme("gleam");  // modern, flat theme that honors background colors

  textbuf = new Fl_Text_Buffer;
  // textbuf->transcoding_warning_action = NULL;
  init_language(argv[0]);
  style_init();
  fl_open_callback(cb);
  fl_message_hotspot(1);

  apply_theme();  // set syntax/background colors for the default (dark) theme

  Fl_Window* window = new_view();

  window->show(1, argv);

  // Re-apply after show(): showing the window re-initializes some scheme colors,
  // so we set the theme colors again here to make sure the dark theme wins.
  apply_theme();

  // #ifndef __APPLE__
  if (argc > 1) {
    load_file(argv[1], -1);
  }
  // #endif

  return Fl::run();
}
