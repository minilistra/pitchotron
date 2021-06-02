# pitchotron
Twibright Pitchotron, a musical transcription program
=====================================================

What does it do?
----------------
Pitchotron takes a media file (audio or video) and extracts the audio track.
Then it analyzes semitones present in music and draws them in a easy to
follow way. The output is a MPEG4 video that contains the original audio track
and the semitone display is scrolling around.

Requirements
------------
* C Compiler
* mplayer with mencoder capable of encoding into MPEG4 (DivX)

Compiling and installation.
---------------------------
Type make. pitchotron-bin and pix-y4m should be compiled. Then change to
root (su -), go back into the same directory and type "make install". Now
pitchotron, pitchotron-bin and pix-y4m are installed in /usr/local/bin. If
/usr/local/bin is not in your $PATH, add it into your $PATH. How this is
done depends on the particular operating system you are using.

Usage
-----
Run "./pitchotron" to see available option and invocation synopsis.

Author, Copyright, Contact
--------------------------
Pitchotron is (c) Karel 'Clock' Kulhavy of Twibright Labs. Pitchotron is a free
software, released under the GNU Public License version 3.0. Send comments
to clock (at) twibright (dot) com.

