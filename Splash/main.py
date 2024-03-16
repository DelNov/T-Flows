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

        self.compile_button = customtkinter.CTkButton(master=self.frame_left, text="Compile Code",
                                              fg_color=("gray75", "gray30"), command=self.compile_code)
        self.compile_button.grid(row=6, column=0, pady=10, padx=20)  # Adjust grid placement as needed

        self.button_2 = customtkinter.CTkButton(master=self.frame_left, text="CTkButton 2",
                                                 fg_color=("gray75", "gray30"), command=self.button_event)
        self.button_2.grid(row=3, column=0, pady=10, padx=20)

        self.button_3 = customtkinter.CTkButton(master=self.frame_left, text="CTkButton 3",
                                                 fg_color=("gray75", "gray30"), command=self.button_event)
        self.button_3.grid(row=4, column=0, pady=10, padx=20)

        self.switch_1 = customtkinter.CTkSwitch(master=self.frame_left)
        self.switch_1.grid(row=9, column=0, pady=10, padx=20, sticky="w")

        self.switch_2 = customtkinter.CTkSwitch(master=self.frame_left, text="Dark Mode",
                                                 command=self.change_mode)
        self.switch_2.grid(row=10, column=0, pady=10, padx=20, sticky="w")

        self.frame_right.rowconfigure((0, 1, 2, 3), weight=1)
        self.frame_right.rowconfigure(7, weight=10)
        self.frame_right.columnconfigure((0, 1), weight=1)
        self.frame_right.columnconfigure(2, weight=0)

        self.frame_info = customtkinter.CTkFrame(master=self.frame_right)
        self.frame_info.grid(row=0, column=0, columnspan=2, rowspan=4, pady=20, padx=20, sticky="nsew")

        self.frame_info.rowconfigure(0, weight=1)
        self.frame_info.columnconfigure(0, weight=1)

        self.label_info_1 = customtkinter.CTkLabel(master=self.frame_info,
                                                    text="CTkLabel: Lorem ipsum dolor sit,\n" +
                                                         "amet consetetur sadipscing elitr,\n" +
                                                         "sed diam nonumy eirmod tempor",
                                                    height=100,
                                                    fg_color=("white", "gray38"),
                                                    justify=tkinter.LEFT)
        self.label_info_1.grid(column=0, row=0, sticky="nwe", padx=15, pady=15)

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

        self.entry = customtkinter.CTkEntry(master=self.frame_right, width=120,
                                            placeholder_text="CTkEntry")
        self.entry.grid(row=8, column=0, columnspan=2, pady=20, padx=20, sticky="we")

        self.button_5 = customtkinter.CTkButton(master=self.frame_right, text="CTkButton",
                                                 command=self.button_event)
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

    def compile_code(self):
        base_path = "../Sources"  # Adjust based on the actual path to the Sources folder relative to the script
        subprograms = ["Process", "Generate", "Divide", "Convert"]
        for subprogram in subprograms:
            try:
                # Construct the path to each subprogram
                subprogram_path = f"{base_path}/{subprogram}"
                # Run the make command in each subprogram's directory
                subprocess.run(["make", "-C", subprogram_path], check=True)
                tkinter.messagebox.showinfo("Compilation Success", f"{subprogram} compiled successfully!")
            except subprocess.CalledProcessError as e:
                tkinter.messagebox.showerror("Compilation Error", f"Failed to compile {subprogram}. Error: {e}")
                break  # Stops the loop if a compilation fails
             
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
