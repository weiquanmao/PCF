#include "melody.h"

const int melody_oneBeat = 400;
const int melody_halfBeat = 200;
const int melody_hold = 128;

void MelodyPlay_Notice()
{
    Beep(note_la, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_do1, melody_oneBeat + melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 3 * melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi, melody_halfBeat);
    Beep(note_mi, melody_halfBeat);

    Beep(note_la, melody_halfBeat + melody_oneBeat);
    Beep(note_so, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
}
void MelodyPlay_CIS()
{
    Beep(note_la, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_do1, melody_oneBeat + melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 3 * melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi, melody_halfBeat);
    Beep(note_mi, melody_halfBeat);

    Beep(note_la, melody_halfBeat + melody_oneBeat);
    Beep(note_so, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_so, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_oneBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold / 2);

    Beep(note_fa, melody_oneBeat + melody_halfBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_fa, melody_halfBeat);
    Beep(note_do1, melody_halfBeat + melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_halfBeat);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold / 2);

    Beep(note_si, melody_halfBeat + melody_oneBeat);
    Beep(note_sfa, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_sfa, melody_oneBeat);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_oneBeat);
    Beep(note_la, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_do1, melody_oneBeat + melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_oneBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(20);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_la, melody_oneBeat + melody_halfBeat);
    Beep(note_so, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_so, 3 * melody_oneBeat);
    Sleep(melody_hold + melody_halfBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold / 2);

    Beep(note_fa, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(20);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_mi1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold + melody_oneBeat);

    Beep(note_do1, melody_oneBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_halfBeat);
    Sleep(20);
    Beep(note_la, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_sso, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_sso, 2 * melody_oneBeat);
    Sleep(melody_hold + melody_oneBeat);
    Beep(note_do1, melody_halfBeat);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_mi1, melody_oneBeat + melody_halfBeat);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_fa1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_re1, 2 * melody_oneBeat);
    Sleep(melody_oneBeat + melody_hold);
    Beep(note_so, melody_halfBeat);
    Sleep(20);
    Beep(note_so, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_do1, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi1, 2 * melody_oneBeat);
    Sleep(melody_hold + 2 * melody_oneBeat);

    Beep(note_la, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_do1, melody_oneBeat + melody_halfBeat);
    Beep(note_so, melody_halfBeat);
    Sleep(20);
    Beep(note_so, melody_oneBeat);
    Sleep(melody_oneBeat + melody_hold);

    Beep(note_fa1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_re1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi1, 4 * melody_oneBeat);

    Beep(note_mi1, melody_oneBeat * 2);
    Sleep(melody_oneBeat + melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_la1, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_so1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_so1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold + melody_halfBeat);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_re1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_halfBeat);
    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_so1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi1, 2 * melody_oneBeat);
    Sleep(melody_hold + melody_oneBeat);
    Beep(note_mi, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_la1, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_so1, 2 * melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi1, melody_halfBeat);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, 2 * melody_oneBeat);
    Sleep(melody_hold + melody_halfBeat);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_re1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_halfBeat);
    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_re1, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_la, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_halfBeat);
    Beep(note_si, melody_halfBeat);

    Beep(note_do1, melody_oneBeat + melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 3 * melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi, melody_halfBeat);
    Beep(note_mi, melody_halfBeat);

    Beep(note_la, melody_halfBeat + melody_oneBeat);
    Beep(note_so, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_so, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_oneBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold / 2);

    Beep(note_fa, melody_oneBeat + melody_halfBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_fa, melody_halfBeat);
    Beep(note_do1, melody_halfBeat + melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_mi, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_halfBeat);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_do1, melody_halfBeat);
    Sleep(melody_hold / 2);

    Beep(note_si, melody_halfBeat + melody_oneBeat);
    Beep(note_sfa, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_sfa, melody_oneBeat);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_oneBeat);
    Beep(note_la, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_do1, melody_oneBeat + melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_mi1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_si, 2 * melody_oneBeat);
    Sleep(melody_hold);
    Sleep(melody_oneBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(20);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold);

    Beep(note_la, melody_oneBeat + melody_halfBeat);
    Beep(note_so, melody_halfBeat);
    Sleep(melody_hold);
    Beep(note_la, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_so, 3 * melody_oneBeat);
    Sleep(melody_hold + melody_halfBeat);
    Beep(note_mi, melody_halfBeat);
    Sleep(melody_hold / 2);

    Beep(note_fa, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_halfBeat);
    Beep(note_si, melody_halfBeat);
    Sleep(20);
    Beep(note_si, melody_oneBeat);
    Sleep(melody_hold);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold);

    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_re1, melody_halfBeat);
    Sleep(20);
    Beep(note_mi1, melody_halfBeat);
    Sleep(melody_hold / 2);
    Beep(note_do1, melody_oneBeat);
    Sleep(melody_hold + melody_oneBeat);

    Beep(note_la, 4 * melody_oneBeat);
}