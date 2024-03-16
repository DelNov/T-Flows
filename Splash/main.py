import tkinter
import tkinter.messagebox
import customtkinter
import sys
import subprocess


customtkinter.set_appearance_mode("Dark")
customtkinter.set_default_color_theme("dark-blue")


class App(customtkinter.CTk):

    WIDTH = 780
    HEIGHT = 520

    def __init__(self):
        super().__init__()

        self.title("SplashFlow complex example")
        self.geometry(f"{App.WIDTH}x{App.HEIGHT}")

        self.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.frame_left = customtkinter.CTkFrame(master=self, width=180, corner_radius=0)
        self.frame_left.grid(row=0, column=0, sticky="nswe")

        self.frame_right = customtkinter.CTkFrame(master=self)
        self.frame_right.grid(row=0, column=1, sticky="nswe", padx=20, pady=20)

        self.frame_left.grid_rowconfigure(0, minsize=10)
        self.frame_left.grid_rowconfigure(5, weight=1)
        self.frame_left.grid_rowconfigure(8, minsize=20)
        self.frame_left.grid_rowconfigure(11, minsize=10)

        self.label_1 = customtkinter.CTkLabel(master=self.frame_left, text="SplashFlows")
        self.label_1.grid(row=1, column=0, pady=10, padx=10)
        self.label_1.configure(font=("Roboto Medium", 16))
    
        self.process_check_var = tkinter.IntVar(value=0)
        self.generate_check_var = tkinter.IntVar(value=0)
        self.divide_check_var = tkinter.IntVar(value=0)
        self.convert_check_var = tkinter.IntVar(value=0)

        self.process_checkbox = customtkinter.CTkCheckBox(master=self.frame_left, text="Process", variable=self.process_check_var)
        self.generate_checkbox = customtkinter.CTkCheckBox(master=self.frame_left, text="Generate", variable=self.generate_check_var)
        self.divide_checkbox = customtkinter.CTkCheckBox(master=self.frame_left, text="Divide", variable=self.divide_check_var)
        self.convert_checkbox = customtkinter.CTkCheckBox(master=self.frame_left, text="Convert", variable=self.convert_check_var)

        self.process_checkbox.grid(row=7, column=0, pady=10, padx=20)
        self.generate_checkbox.grid(row=8, column=0, pady=10, padx=20)
        self.divide_checkbox.grid(row=9, column=0, pady=10, padx=20)
        self.convert_checkbox.grid(row=10, column=0, pady=10, padx=20)

    
##        self.compile_button = customtkinter.CTkButton(master=self.frame_left, text="Compile Code",
##                                              fg_color=("gray75", "gray30"), command=self.compile_code)
##        self.compile_button.grid(row=6, column=0, pady=10, padx=20)  # Adjust grid placement as needed
        
        self.compile_button = customtkinter.CTkButton(master=self.frame_left, text="Compile Code",
                                                 command=self.compile_code)
        self.compile_button.grid(row=6, column=0, pady=10, padx=20)
        
        # Switch for changing program theme
        self.switch_2 = customtkinter.CTkSwitch(master=self.frame_left, text="Dark Mode",
                                                 command=self.change_mode)
        self.switch_2.grid(row=12, column=0, pady=10, padx=20, sticky="w")

        self.frame_right.rowconfigure((0, 1, 2, 3), weight=1)
        self.frame_right.rowconfigure(7, weight=10)
        self.frame_right.columnconfigure((0, 1), weight=1)
        self.frame_right.columnconfigure(2, weight=0)

        self.frame_info = customtkinter.CTkFrame(master=self.frame_right)
        self.frame_info.grid(row=0, column=0, columnspan=2, rowspan=4, pady=20, padx=20, sticky="nsew")

        self.frame_info.rowconfigure(0, weight=1)
        self.frame_info.columnconfigure(0, weight=1)
        
        # Create a Text widget for displaying the compilation output (text box)
        self.output_text = tkinter.Text(master=self.frame_info, height=20, width=80, bg="black", fg="white", wrap=tkinter.WORD)
        self.output_text.grid(row=0, column=0, sticky="nsew", padx=15, pady=15)

        # Create a Scrollbar and attach it to the Text widget
        self.output_scrollbar = tkinter.Scrollbar(master=self.frame_info, command=self.output_text.yview)
        self.output_scrollbar.grid(row=0, column=1, sticky='nsew')
        self.output_text.config(yscrollcommand=self.output_scrollbar.set)

        # Adjust the frame_info grid configuration to accommodate the scrollbar
        self.frame_info.grid_columnconfigure(0, weight=1)
        self.frame_info.grid_columnconfigure(1, weight=0)
        self.frame_info.grid_rowconfigure(0, weight=1)

        self.progressbar = customtkinter.CTkProgressBar(master=self.frame_info)
        self.progressbar.grid(row=1, column=0, sticky="ew", padx=15, pady=15)

        self.radio_var = tkinter.IntVar(value=0)

        self.label_radio_group = customtkinter.CTkLabel(master=self.frame_right,
                                                        text="CTkRadioButton Group:")
        self.label_radio_group.grid(row=0, column=2, columnspan=1, pady=20, padx=10, sticky="")

        self.radio_button_1 = customtkinter.CTkRadioButton(master=self.frame_right,
                                                           variable=self.radio_var,
                                                           value=0)
        self.radio_button_1.grid(row=1, column=2, pady=10, padx=20, sticky="n")

        self.radio_button_2 = customtkinter.CTkRadioButton(master=self.frame_right,
                                                           variable=self.radio_var,
                                                           value=1)
        self.radio_button_2.grid(row=2, column=2, pady=10, padx=20, sticky="n")

        self.radio_button_3 = customtkinter.CTkRadioButton(master=self.frame_right,
                                                           variable=self.radio_var,
                                                           value=2)
        self.radio_button_3.grid(row=3, column=2, pady=10, padx=20, sticky="n")

        self.slider_1 = customtkinter.CTkSlider(master=self.frame_right, from_=0, to=1,
                                                 number_of_steps=3, command=self.progressbar.set)
        self.slider_1.grid(row=4, column=0, columnspan=2, pady=10, padx=20, sticky="we")

        self.slider_2 = customtkinter.CTkSlider(master=self.frame_right, command=self.progressbar.set)
        self.slider_2.grid(row=5, column=0, columnspan=2, pady=10, padx=20, sticky="we")

        self.slider_button_1 = customtkinter.CTkButton(master=self.frame_right, height=25,
                                                        text="CTkButton", command=self.button_event)
        self.slider_button_1.grid(row=4, column=2, columnspan=1, pady=10, padx=20, sticky="we")

        self.slider_button_2 = customtkinter.CTkButton(master=self.frame_right, height=25,
                                                        text="CTkButton", command=self.button_event)
        self.slider_button_2.grid(row=5, column=2, columnspan=1, pady=10, padx=20, sticky="we")

        self.checkbox_button_1 = customtkinter.CTkButton(master=self.frame_right, height=25,
                                                          text="CTkButton", border_width=3,
                                                          fg_color=None, command=self.button_event)
        self.checkbox_button_1.grid(row=6, column=2, columnspan=1, pady=10, padx=20, sticky="we")

        self.check_box_1 = customtkinter.CTkCheckBox(master=self.frame_right, text="CTkCheckBox")
        self.check_box_1.grid(row=6, column=0, pady=10, padx=20, sticky="w")

        self.check_box_2 = customtkinter.CTkCheckBox(master=self.frame_right, text="CTkCheckBox")
        self.check_box_2.grid(row=6, column=1, pady=10, padx=20, sticky="w")


        # Plan B: save me with the terminal :/
        self.entry = customtkinter.CTkEntry(master=self.frame_right, width=120,
                                    placeholder_text="top")
        self.entry.grid(row=8, column=0, columnspan=2, pady=20, padx=20, sticky="we")

        self.button_5 = customtkinter.CTkButton(master=self.frame_right, text="Execute",
                                                command=self.execute_command)
        self.button_5.grid(row=8, column=2, columnspan=1, pady=20, padx=20, sticky="we")



        self.radio_button_1.select()
        self.switch_2.select()
        self.slider_1.set(0.2)
        self.slider_2.set(0.7)
        self.progressbar.set(0.5)
        self.slider_button_1.configure(state=tkinter.DISABLED, text="Disabled Button")
        self.radio_button_3.configure(state=tkinter.DISABLED)
        self.check_box_1.configure(state=tkinter.DISABLED, text="CheckBox disabled")
        self.check_box_2.select()

    def execute_command(self):
        # Retrieve the command from the entry widget
        command = self.entry.get()
        if not command:  # Fallback if the entry is empty, use the placeholder or a default command
            command = "top"
        
        # Platform-specific logic to open a terminal and run the command
        if sys.platform.startswith('linux'):
            # For Linux, using gnome-terminal as an example
            subprocess.Popen(['gnome-terminal', '--', 'bash', '-c', command])
        elif sys.platform.startswith('win32'):
            # For Windows
            subprocess.Popen(['cmd', '/c', command])
        elif sys.platform.startswith('darwin'):
            # For macOS
            subprocess.Popen(['open', '-a', 'Terminal', command])
        else:
            print(f"Unsupported platform: {sys.platform}")
            
        
    def compile_code(self):
        self.output_text.delete('1.0', tkinter.END)  # Clear the Text widget at the start

        subprograms = {
            "Process": self.process_check_var,
            "Generate": self.generate_check_var,
            "Divide": self.divide_check_var,
            "Convert": self.convert_check_var,
        }
        selected_subprograms = {name: var.get() for name, var in subprograms.items() if var.get() == 1}

        if not any(selected_subprograms.values()):
            tkinter.messagebox.showerror("Error", "No choice was detected, please choose the part of the code you want to compile.")
            return

        compile_success = True  # Flag to track overall compilation success
        base_path = "../Sources"
        for subprogram, selected in selected_subprograms.items():
            if selected:
                subprogram_path = f"{base_path}/{subprogram}"
                try:
                    process = subprocess.Popen(["make", "-C", subprogram_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                    for line in iter(process.stdout.readline, ''):
                        self.output_text.insert(tkinter.END, line)
                        self.output_text.see(tkinter.END)
                        self.update_idletasks()
                    process.wait()
                    if process.returncode != 0:
                        tkinter.messagebox.showerror("Compilation Error", f"Failed to compile {subprogram}.")
                        compile_success = False
                        break
                except Exception as e:
                    tkinter.messagebox.showerror("Compilation Error", f"Failed to compile {subprogram}. Error: {e}")
                    compile_success = False
                    break

        if compile_success:
            self.output_text.insert(tkinter.END, "Compilation successful for all selected subprograms!")
        
             
    def button_event(self):
        print("Button pressed")

    def change_mode(self):
        if self.switch_2.get() == 1:
            customtkinter.set_appearance_mode("dark")
        else:
            customtkinter.set_appearance_mode("light")

    def on_closing(self, event=0):
        self.destroy()

    def start(self):
        self.mainloop()


if __name__ == "__main__":
    app = App()
    app.start()
