import os, sys

if __name__ == '__main__':
	main()

	def installMac(shell, name, download, folder, website, origName, extract=True, build="", makefile=""):
		print "Checking", name, "..."
		try:
			subprocess.call([shell])
		except:
			print "...you don't have", name, "MUSCLE installed..."
			fileDwn = urllib.urlretrieve(download)
			if extract:
				subprocess.call(["tar","-zxf", fileDwn[0]])
			if build:
				subprocess.call(['make', "-C " + build, "-f" + makefile])
				subprocess.call(['mv', pathName + "raxmlHPC-SSE3", "raxml"])
		
			subprocess.call(["mv", folder + "/" + origName, os.getcwd() + "/" + shell])
			print "...testing", name, "..."
			try:
				subprocess.call([shell])
			except:
				print "!!!Can't install", name, ". Sorry! Loading up its website, and exiting."
				webrowser.open(website)
				sys.exit()
			print "...", name, "installed!"

def main():
	print "\n\nWelcome to phyloGenerator's installer!"
	print "---Please go to http://willpearse.github.com/phyloGenerator for help"
	print "---Written by Will Pearse (will.pearse@gmail.com)"
	
	if sys.platform == 'darwin':
		if os.getuid() is not 0:
			print "\nSorry! This script needs admin privileges! Please run again like this:"
			print "'sudo python installer.py'"
			print "- then key in your password"
			print "...exiting..."
			sys.exit()
		
		#Python
		print "Checking python..."
		if sys.version_info.minor >= 6 and sys.version_info.major==2:
			print "...Python version is OK!"
		else:
			print "...your version of Python isn't recent enough, and I can't change that from within Python..."
			print "\nI'm opening a browser window to install a newer version of Python. If that doesn't work, go to http://www.python.org/getit/ and install Pytohn 2."
			webrowser.open("http://www.python.org/getit/")
			sys.exit()
		
		#GCC
		print "Checking C/C++ compilers..."
		try:
			subprocess.call(['make'])
		except:
			print "...you don't have a C compiler installed. I can't fix this from Python - on a Mac, you need to install 'Xcode'"
			print "\nI'm opening a browser window to help install it. If that doesn't work, go to http://developer.apple.com/xcode/ and install it"
			webrowser.open("http://developer.apple.com/xcode/")
			sys.exit()
		print "...compilers are installed!"
	
		#NumPy
		print "Checking for Numerical Python (NumPy)..."
		try:
			import numpy
			"...NumPy already installed!"
		except:
			print "...you don't have NumPy installed. I'm going to download it for you as a disk image. Please just open it up once it's downloaded, and double-click the installer."
		webbrowser.open("http://downloads.sourceforge.net/project/numpy/NumPy/1.6.1/numpy-1.6.1-py2.6-python.org-macosx10.3.dmg?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fnumpy%2Ffiles%2F&ts=1326294925")
			print "I'm going to exit now. Re-start this once you've installed NumPy"
			sys.exit()
		
		#SciPy
		print "Checking for Scientific Python (SciPy)..."
		try:
			import scipy
			"...SciPy already installed!"
		except:
			print "...you don't have SciPy installed. I'm going to download it for you as a disk image. Please just open it up once it's downloaded, and double-click the installer."
		webbrowser.open("http://downloads.sourceforge.net/project/scipy/scipy/0.10.0/scipy-0.10.0-py2.6-python.org-macosx10.3.dmg?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fscipy%2Ffiles%2F&ts=1326295454")
			print "I'm going to exit now. Re-start this once you've installed SciPy"
			sys.exit()
		
		#BioPython
		print "Checking for BioPython..."
		try:
			from Bio import SeqIO
			"...BioPython already installed!"
		except:
			fileDwn = urllib.urlretrieve(http://biopython.org/DIST/biopython-1.58.tar.gz)
			subprocess.call(["tar","-zxf", fileDwn[0]])
			filDwn = list(fileDwn[0].partition("/"))
			path = fileDwn.pop()
			subprocess.call(["python", "biopython-1.58/setup.py", "build"])
			sudo subprocess.call(["sudo", "python", "install"])
			print "I'm going to exit now. BioPython should have been successfully installed, but I can't check until I reload Python."
			print " - just restart this script, and all should be fine."
			sys.exit()
		
		#MUSCLE
		installMac('muscle', "MUSCLE", "http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86darwin32.tar.gz", "", "http://www.drive5.com/muscle/", "muscle3.8.31_i86darwin32", extract=True, build=False)
		
		#MAFFT
		print "Checking MAFFT..."
		try:
			subprocess.call(['mafft'])
		except:
			print "...you don't have MAFFT installed..."
			fileDwn = urllib.urlretrieve("http://mafft.cbrc.jp/alignment/software/mafft-6.864-mac-universal.mpkg.zip")
			subprocess.call(["mv", fileDwn[0], os.getcwd()+"/doubleClickMe"])
			print "\nIn your current working directory there's a file called 'doubleClickMe'. Double click it to launch the mafft installer. When you're done, just hit enter"
			waitInput = raw_input("Installing MAFFT...")
			print "...testing MAFFT..."
			try:
				subprocess.call(['mafft'])
			except:
				print "!!!Something's gone wrong with your MAFFT install. Sorry! Loading up its website, and exiting."
				webrowser.open("http://mafft.cbrc.jp/alignment/software/")
				sys.exit()
			print "...MAFFT installed!"
		
		#Clustal-O
		installMac('clustalo', "Clustal Omega", "http://www.clustal.org/omega/clustal-omega-1.0.3-Mac-ppc-x86_64", "", "http://www.clustal.org/omega/", "clustalo", extract=True, build=False)
	
		#RAxML
		print "Checking RAxML..."
		try:
			subprocess.call(['raxml'])
		except:
			print "...you don't have RAxML installed..."
			print "Downloading RAxML (this will use your web browser, sorry)"
			webbrowser.open("https://github.com/stamatak/standard-RAxML/zipball/master")
			print "...done.\nPlease type in the *directory and name* of the zip file your browser has just downloaded. For example, '/Users/will/Downloads/stamatak-standard-RAxML-f8cdf42.zip'."
			print " - remember, you can just drag and drop the file into Terminal and your computer will fill in the address!"
			locker = True
			while locker:
				zipFile = raw_input("")
				try:
					subprocess.call(['zip', '-xvf', zipFile, "RAxML"])
					locker = False
				except:
					print "Couldn't load file", zipFile, " - please try again"
			zipFile = list(zipFile.partition("/"))
			fileName = zipFile.pop()
			pathName = "".join(fileName)
			subprocess.call(['make', "-C" + pathName, '-f Makefile.SSE3.gcc'])
			subprocess.call(['mv', pathName + "raxmlHPC-SSE3", "raxml"])
			print "...testing RAxML..."
			try:
				subprocess.call(['raxml'])
			except:
				print "!!!Something's gone wrong with your RAxML (SSE3) install. Sorry! Loading up its website, and exiting."
				webrowser.open("https://github.com/stamatak/standard-RAxML")
				sys.exit()
			print "...RAxML (SSE3) installed!"
		
		#trimAl
		installMac('trimal', "trimAl", "http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz", "", "http://trimal.cgenomics.org/", "trimAl", extract=True, build="source", makefile="makefile")
		
		#PATHd8
		installMac('PATHd8', 'PATHd8', 'http://www2.math.su.se/PATHd8/mac/PATHd8', 'http://www2.math.su.se/PATHd8/', 'PATHd8')
		
		print "Congratulations! You've installed *everything*!"
		print "...now just run phyloGenerator.py ('python phyloGenerator.py') and get going!"
	
	else:
		print "Sorry, but this installer only works on a Mac at the moment. Good luck! :p"
		sys.exit()
