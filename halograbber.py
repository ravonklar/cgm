import boxsdk as box
import sys
import os

oauth = box.OAuth2(
    client_id=sys.argv[1],
    client_secret=sys.argv[2]
)

auth_url, csrf_token = oauth.get_authorization_url(None)

print(auth_url)
#I HAVE NO IDEA HOW THE CODES GET IN HERE. MAY HAVE TO COPY AND PASTE.
print('WAITING FOR AUTH CODE')
auth_code = input()
#print('WAITING FOR TOKEN')
#csrf_token_1 = input()
#assert csrf_token_1 == csrf_token

access_token, refresh_token = oauth.authenticate(auth_code)
client = box.Client(oauth)

user = client.user().get()
print(f'User ID is {user.id}')

folder_id = 252615424304
#file_list = os.listdir('/scratch/ravonklar/')
items = client.folder(folder_id).get_items()
for halofile in items:
    name = halofile.name
    with open(f'TNG100-1/cutouts/{name}', 'wb') as open_file:
        halofile.download_to(open_file)
        print(f'Uploaded {name}')
        open_file.close()
